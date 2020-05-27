# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 15:36:06 2020

@author: LucasTsubasaYang
"""
import re, sys

sam_file = sys.argv[1]
ref_path = sys.argv[2]
kmer = int(sys.argv[3])
filtered_sam_path = sys.argv[4]

seq = ''
f = open(ref_path)
for line in f:
    if not line.startswith('>'):
        seq += line.strip()

comp = re.compile('[0-9]+[SHMNDIPX]')
path = sam_file
f = open(path)
g = open(filtered_sam_path, 'w')
# as for hisat2 sam file only delete 'N' in CIGAR
# and change the startpos by adding N span length
for line in f:
    if line.startswith('@'):
        continue
    cols = line.strip().split('\t')
    start_pos = int(cols[3])
    current_pos = start_pos
    read_current_pos = 1
    for item in re.findall(comp, cols[5]):
        span, letter = item[:-1], item[-1]
        if letter in 'MND':
            current_pos += int(span)
        if letter in 'MIS':
            read_current_pos += int(span)
        if letter == 'N' and int(span) > 20000 and 10 < start_pos < 150:
            gap_len = int(span)
            read_current_index = read_current_pos-1
            current_index = current_pos-1
            near_query_seq = cols[9][read_current_index-min(30, read_current_index):read_current_index+kmer+30]
            near_genome_region = seq[current_index-30:current_index+kmer+30]
            if near_query_seq:
                print('\t'.join([cols[0], cols[1], cols[3], str(current_pos), near_query_seq, near_genome_region]))
                new_cols = cols
                new_cols[3] = str(start_pos + gap_len)
                new_cols[5] = ''.join([i for i in re.findall(comp, cols[5]) if 'N' not in i])
                g.write('\t'.join(new_cols[:11]+[str(current_pos)])+'\n')
            break
f.close()
g.close()