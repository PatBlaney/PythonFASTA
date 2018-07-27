#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 22:20:00 2018



@author: pat
"""

# Open the fasta file for analysis or return error message if the file can't 
# be opened
try:
    fasta_file = open('dna2.fasta')
except IOError:
    print('File does not exist')

# Loop through each line in the fasta file and count the total number of 
# records
record_count = 0 
for entry in fasta_file:
    
    # Seach each line for '>' character that always begins a new DNA
    # sequence record 
    if entry[0] == '>':
        
        # Add to the total record count for each found ">" character
        record_count += 1

# Print out the record count
print('The total number of records is %d' % (record_count))

# Use the DNA sequence headers to identify the various lengths of the sequences
# Start by creating dictionary that connects the length of the DNA sequence
# with the sequence itself
# Return to the first line of the file for reading again
fasta_file.seek(0)

# Initalize the empty dictionary
seq_info = {}

for line in fasta_file:
    
    # Remove end of line characters for easier string parsing
    line = line.rstrip()

    # Check if the line of the file is a header
    if line[0] == '>':
        
        # Split the header into individual strings and isolate the sequence
        # length and use this as the key in the dictionary
        header = line.split()
        numeric_metadata = header[0].split('|')
        seq_identifier = numeric_metadata[4]
        seq_info[seq_identifier] = ''
    
    # Store the sequence under its respective sequence length
    else:
        seq_info[seq_identifier] = seq_info[seq_identifier] + line

# Print out all sequence lengths
print('The length of the DNA sequences in this file are...')

max_seq_length = 0
min_seq_length = 1000000
for seq_string in seq_info.values():

    print(len(seq_string))

    if len(seq_string) > max_seq_length:
        max_seq_length = len(seq_string)

    if len(seq_string) < min_seq_length:
        min_seq_length = len(seq_string)

print('The longest DNA sequence in the file is %d' % (max_seq_length))
print('The shortest DNA sequence in the file is %d' % (min_seq_length))

# Use the dictionary of DNA sequences to identify all ORFs within each sequence
# Ask user to define valid reading frame 
reading_frame = input('Please input the desired reading frame. Select either 1, 2, or 3. ')
print('Current reading frame --> %s' % (reading_frame))

# Create function that will find the ORFs for any of the 3 reading frames
def open_reading_frames(fasta_seq, frame = 1):
    """Determies all ORFs of a given sequence for the specified reading frame,
    then returns all ORFs as a list of strings"""
    start_codon = 'ATG'
    stop_codons = ['TAG', 'TGA', 'TAA']
    
    orf_start_positions = []
    orf_stop_positions = []
    
    for i in range((int(frame) - 1), len(fasta_seq), 3):
        codon = fasta_seq[i:i + 3].upper()
        
        if codon == start_codon:
            orf_start_positions.append(i + 1)
        if codon in stop_codons:
            orf_stop_positions.append(i + 1)
    
    if len(orf_start_positions) == 0 or len(orf_stop_positions) == 0:
        print('This DNA sequence has no ORFs for this reading frame')
    else:
        print('ORF Start Codon Index --> %d' % (orf_start_positions[0]))
        print('ORF Stop Codon Index --> %d' % (orf_stop_positions[-1]))
        orf_length = orf_stop_positions[-1] - orf_start_positions[0]
        if orf_length > 0:
            print('The longest possible ORF for this DNA sequence in this reading frame is %d nucleotides' % (orf_length))
        else:
            print('This DNA sequence has no ORFs for this reading frame')
        
# Loop through all the values of the DNA sequence dictionary
for dna_seq in seq_info.values():
    open_reading_frames(dna_seq, reading_frame)
    
# Allow user to input string of target repeat pattern
repeat_pattern = input('Please input the target repeat pattern. ').upper()
print('Current target repeat pattern --> %s' % (repeat_pattern))

def identify_repeats(seq, pattern):
    """Find and count all repeats of a given pattern in a DNA sequence"""
    count_of_repeats = 0
    for i in range(0, len(dna_seq), 1):
        dna_slice = seq[i:(i + len(pattern))]
        
        if dna_slice == pattern:
            count_of_repeats += 1
            
    print('The sequence %s is repeated %d times for this DNA sequence' % (pattern, count_of_repeats))

for dna_seq in seq_info.values():
    identify_repeats(dna_seq, repeat_pattern)