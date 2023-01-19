#! /usr/bin/env python3

"""
@author: Lawrence Gordon
"""

import argparse
import csv
from copy import deepcopy
import re
from Bio.Seq import Seq
#import Bio.Data.IUPACData as bdi
from itertools import product
import random
import numpy as np
import textwrap

# parses command-line arguments
def parse_args():
    parser = argparse.ArgumentParser(description="Synthetic Control Design")
    parser.add_argument("manifest", help="csv formatted manifest with sequences and positions")
    parser.add_argument("--seq_length", default=500, type=int, help="length of the desired sequence")
    parser.add_argument("--gc", default=0.50, type=float, help="gc content of the desired sequence")
    parser.add_argument("--ann_length", default=6, type=int, help="length of the terminal primer sequence to consider")
    parser.add_argument("--h_length", default=5, type=int, help="length of the homopolymers to avoid; zero or one will allow any homopolymers")
    parser.add_argument("--out", "-o", default="synthetic_control.fasta", help="output filename")
    return parser.parse_args()

# parse manifest
def parse_manifest(manifest):
    # [{'Name': '16S', ... 'reverse-start-pos': '431'}, ...]
    manifest_data = []

    try:
        with open(manifest) as manifest_file:
            reader = csv.DictReader(manifest_file)
            for line in reader:
                manifest_data.append(line)
    except IOError:
        print("Unable to open file {0}".format(manifest))

    return manifest_data

# perform reverse complementing of manifest data:
def reverse_complement(manifest_data):
    complemented_manifest_data = deepcopy(manifest_data)

    for primer in complemented_manifest_data:
        # change I as K prior to reverse complement
        reverse_primer = Seq(re.sub("I", "K", primer['Reverse']))
        reverse_revcomp = str(reverse_primer.reverse_complement())
        primer['Reverse'] = reverse_revcomp
    
    return complemented_manifest_data

# disambiguate sequences
def disamb_sequences(seq_list):
    """
    ambiguous_dna_complement = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
    "M": "K",
    "R": "Y",
    "W": "W",
    "S": "S",
    "Y": "R",
    "I": "M",
    "K": "M",
    "V": "B",
    "H": "D",
    "D": "H",
    "B": "V",
    "X": "X",
    "N": "N",
    }
    """

    ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "I": "GT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
    }

    #d = bdi.ambiguous_dna_values
    disamb_seq_list = []

    # disambiguate all sequences
    for seq in seq_list:
        for disamb_seq in list(map("".join, product(*map(ambiguous_dna_values.get, seq)))):
            disamb_seq_list.append(disamb_seq)

    return disamb_seq_list

# set parameters for randomized sequence
def derive_parameters(complemented_manifest_data, annealing_length):
    # define set of invalid string
    invalid_seqs = set()

    for primer in complemented_manifest_data:
        if len(primer['Reverse']) != 0:
            # take 3' end of forward primer
            forward_ann = primer['Forward'][-annealing_length:]
            invalid_seqs.add(forward_ann)
            # take 5' end of revcomp reverse primer
            reverse_ann = primer['Reverse'][:annealing_length]
            invalid_seqs.add(reverse_ann)
        else: # condition to add all of rxn enzyme
            invalid_seqs.add(primer['Forward'][-annealing_length:]) ### new addition for deterministic
            ### trim rxn enzyme to annealing length for consistency

    return set(disamb_sequences(invalid_seqs))

# add desired length homopolymers to prevent from being added to the sequence
def invalid_homopolymers(invalid_seqs, homopolymer_length):
    invalid_seqs = deepcopy(invalid_seqs)
    homopolymers = set([nucleotide*homopolymer_length for nucleotide in "ACGT"])
    return invalid_seqs.union(homopolymers)    

# some sequences might overlap in primer regions
def get_valid_sequences(complemented_manifest_data):
    valid_seqs = ""

    for primer in complemented_manifest_data:
        valid_seqs += primer['Forward']
        valid_seqs += primer['Reverse']

    return valid_seqs

# generate dictionary of keys --> valid ambiguous sequences
# values --> list of unambiguous valid sequences
def get_valid_sequences2(complemented_manifest_data):
    valid_seqs = {primer[direction]:disamb_sequences([primer[direction]]) \
        for primer in complemented_manifest_data \
        for direction in ("Forward", "Reverse")}

    return valid_seqs

# function for overwriting string with sequence
def insert_sequence(original_sequence, insert, start):
    end = start + len(insert)
    new_sequence = original_sequence[:start] + \
                        insert + \
                        original_sequence[end:]
    
    return new_sequence

# deterministic sequence design
def deterministic_sequence_design(complemented_manifest_data, desired_gc, desired_length, invalid_seqs, valid_seqs, annealing_length):
    nucleotides = "GCAT"
    nucleotides_list = ["G", "C", "A", "T"]
    nucleotides_set = set(["G", "C", "A", "T"])

    # define chance to pick nucleotide
    gc_chance = desired_gc/2
    at_chance = (1 - desired_gc)/2
    nucleotide_weights = {"G": gc_chance, "C": gc_chance, "A": at_chance, "T": at_chance}
    weights = [gc_chance, gc_chance, at_chance, at_chance]
    sequence_pass = False
    
    # generate dictionary of keys --> invalid substrings[:-1] and 
    # values --> (list of valid nuc, weights)
    valid_continuation = {}
    for seq in invalid_seqs:
        sub = seq[:-1]
        last = set(seq[-1])
        valid_last = list(nucleotides_set - last)
        valid_weights = [nucleotide_weights[nuc] for nuc in valid_last]
        valid_continuation[sub] = (valid_last, valid_weights)

    # populate control sequence with desired sequences
    control_sequence = "X" * desired_length
    for primer in complemented_manifest_data:
            # insert forward primer
            forward_seq = primer['Forward']
            f_pos = int(primer['forward-start-pos'])-1
            control_sequence = insert_sequence(control_sequence, forward_seq, f_pos)

            # insert reverse primer
            reverse_seq = primer['Reverse']
            if reverse_seq != "":
                r_pos = int(primer['reverse-start-pos'])-1
                control_sequence = insert_sequence(control_sequence, reverse_seq, r_pos)
    
    starter_sequence = deepcopy([*control_sequence])


    print(invalid_seqs)
    print(valid_seqs)
    print(valid_continuation)

    iterations = 0
    while not sequence_pass:
        iterations += 1; print("Iteration {0}".format(iterations))
        control_sequence = deepcopy(starter_sequence)

        # main loop for filling in rest of control sequence
        for nucleotide in range(desired_length):
            if control_sequence[nucleotide] == "X":
                prev_seq = ("").join(control_sequence[nucleotide-annealing_length+1:nucleotide])
                valid_continuation.setdefault(prev_seq, (nucleotides, weights))
                curr_nucs, curr_weights = valid_continuation[prev_seq]
                control_sequence[nucleotide] = random.choices(curr_nucs, curr_weights, k=1)[0]

        control_sequence = ("").join(control_sequence)

        # begin searching for invalid sequences
        sequence_pass = True
        for seq in invalid_seqs:
            control_count = control_sequence.count(seq)
            valid_count = 0; appearances = 0
            for valid_amb, valid_disambs in valid_seqs.items():
                appearances = max([valid_disamb.count(seq) for valid_disamb in valid_disambs])
                valid_count += appearances
                #max([valid_disamb.count(seq) for valid_disamb in valid_disambs])

                # need to check if counter increased from ambiguous sequence
                if appearances != 0 and valid_amb.count(seq) == 0:
                    control_count += 1

                appearances = 0
                print(valid_amb, seq, control_count, valid_count)
            #print(seq, control_sequence.count(seq), valid_count)
            if control_count > valid_count:
                sequence_pass = False
                break

    return control_sequence

    """
    iterations = 0
    while not sequence_pass:
        # generate random sequence
        control_sequence = ("").join(random.choices(nucleotides, weights=weights, k=desired_length))
        iterations+=1
        
        # populate sequence with manifest data
        for primer in complemented_manifest_data:
            # insert forward primer
            forward_seq = primer['Forward']
            f_pos = int(primer['forward-start-pos'])-1
            control_sequence = insert_sequence(control_sequence, forward_seq, f_pos)

            # insert reverse primer
            reverse_seq = primer['Reverse']
            if reverse_seq != "":
                r_pos = int(primer['reverse-start-pos'])-1
                control_sequence = insert_sequence(control_sequence, reverse_seq, r_pos)

        #print(control_sequence) if iterations % 100 == 0 else None

        # begin searching for invalid sequences
        sequence_pass = True
        for seq in invalid_seqs:
            if control_sequence.count(seq) > valid_seqs.count(seq):
                sequence_pass = False
                
        # future: perform a check for tolerable differences in GC content

    return control_sequence
    """

# generate sequence with desired features
def sequence_design(complemented_manifest_data, desired_gc, desired_length, invalid_seqs, valid_seqs):
    # define chance to pick nucleotide
    nucleotides = "GCAT"
    gc_chance = desired_gc/2
    at_chance = (1 - desired_gc)/2
    weights = [gc_chance, gc_chance, at_chance, at_chance]
    sequence_pass = False

    iterations = 0
    while not sequence_pass:
        # generate random sequence
        control_sequence = ("").join(random.choices(nucleotides, weights=weights, k=desired_length))
        iterations+=1
        
        # populate sequence with manifest data
        for primer in complemented_manifest_data:
            # insert forward primer
            forward_seq = primer['Forward']
            f_pos = int(primer['forward-start-pos'])-1
            control_sequence = insert_sequence(control_sequence, forward_seq, f_pos)

            # insert reverse primer
            reverse_seq = primer['Reverse']
            if reverse_seq != "":
                r_pos = int(primer['reverse-start-pos'])-1
                control_sequence = insert_sequence(control_sequence, reverse_seq, r_pos)

        #print(control_sequence) if iterations % 100 == 0 else None

        # begin searching for invalid sequences
        sequence_pass = True
        for seq in invalid_seqs:
            if control_sequence.count(seq) > valid_seqs.count(seq):
                sequence_pass = False
                
        # future: perform a check for tolerable differences in GC content

    return control_sequence

# get standard nucleotide compsition for a sequence
def nucleotide_composition(sequence):
    seq_len = len(sequence)
    G = round(((sequence.count("G")/seq_len)*100), 2)
    C = round(((sequence.count("C")/seq_len)*100), 2)
    A = round(((sequence.count("A")/seq_len)*100), 2)
    T = round(((sequence.count("T")/seq_len)*100), 2)

    return (G, C, A, T)

# print sequences with primers highlighted
def print_sequences(control_sequence, complemented_manifest_data, header, out_file):

    html_colors = ["blue", "fuchsia", "green", "maroon", "red", "teal", "yellow"]
    num_colors = len(html_colors)
    
    sequence_list = [*control_sequence]

    # loop to color highlight primer sequences
    for idx, primer in enumerate(complemented_manifest_data):
        color = html_colors[idx%num_colors]

        # first color forward data
        start_idx = int(primer['forward-start-pos']) - 1
        end_idx = start_idx + len(primer['Forward']) - 1
        # add proper html format to color text
        sequence_list[start_idx] = '<span style="color:' + color + '">' + sequence_list[start_idx]
        sequence_list[end_idx] = sequence_list[end_idx] + "</span>"

        # then match reverse
        if len(primer['Reverse']) != 0:
            start_idx = int(primer['reverse-start-pos']) - 1
            end_idx = start_idx + len(primer['Reverse']) - 1
            sequence_list[start_idx] = '<span style="color:' + color + '">' + sequence_list[start_idx]
            sequence_list[end_idx] = sequence_list[end_idx] + "</span>"


    # loop to format string as fasta file
    for idx in range(len(sequence_list)):
        if (idx + 1) % 50 == 0:
            sequence_list[idx] += f" ---{idx+1}</br>"
        elif (idx + 1) % 10 == 0:
            sequence_list[idx] += " "

    # name out_file with html syntax
    out_file = (".").join(out_file.split(".")[:-1])+".html" if not out_file.endswith(".html") else out_file

    # write out html file
    with open(out_file, "w") as html:
        html.write("<tt>" + header + "</br>")
        html.write(("").join(sequence_list))
        html.write("<tt>")

# final sequence write out as fasta file
def write_sequence(control_sequence, out_file):
    # generate nucleotide counts
    G, C, A, T = nucleotide_composition(control_sequence)

    # text wrapping for fasta
    wrapper = textwrap.TextWrapper(width=50)
    seq = wrapper.fill(text=control_sequence)

    header = ">synthetic primer control len={0} G:{1}% C:{2}% A:{3}% T:{4}%".format(len(control_sequence), G, C, A, T,seq)

    with open(out_file, "w") as out:
        out.write("{0}\n{1}\n".format(header, seq))

    return header

def main():
    args = parse_args()
    manifest_data = parse_manifest(args.manifest)
    complemented_manifest_data = reverse_complement(manifest_data)
    invalid_seqs = derive_parameters(complemented_manifest_data, args.ann_length)

    # add invalid homopolymers
    if args.h_length != 0 or args.h_length != 1:
        invalid_seqs = invalid_homopolymers(invalid_seqs, args.h_length)
    #valid_seqs = get_valid_sequences(complemented_manifest_data)
    valid_seqs2 = get_valid_sequences2(complemented_manifest_data)

    test = []
    for lst in valid_seqs2.values():
        for val in lst:
            test.append(val)
    
    print(f"Valid sequences: {len(valid_seqs2)}\n{valid_seqs2}")
    ###control_sequence = sequence_design(complemented_manifest_data, args.gc, args.seq_length, invalid_seqs, valid_seqs)
    control_sequence = deterministic_sequence_design(complemented_manifest_data, args.gc, args.seq_length, invalid_seqs, valid_seqs2, args.ann_length)
    header = write_sequence(control_sequence, args.out)
    print_sequences(control_sequence, complemented_manifest_data, header, args.out)

main()
