# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:06:38 2020

@author: LucasTsubasaYang
"""

import sys, re

sam_path = sys.argv[1]
ref_path = sys.argv[2]
kmer = int(sys.argv[3])
filtered_sam_path = sys.argv[4]
spliced_read_path = sys.argv[5]

seq = ''
f = open(ref_path)
for line in f:
    if not line.startswith('>'):
        seq += line.strip()

# find subgenome pos by split reads 
# shortage: it cannot find the last subgenome
comp = re.compile('[0-9]*[SHMNDIPX]')
f = open(sam_path)
g = open(filtered_sam_path, 'w')
h = open(spliced_read_path, 'w')
pre_5_items = []
post_5_items = []
for line in f:
    if line.startswith('@'):
        continue
    cols = line.strip().split('\t')
    start_pos = int(cols[3])
    current_pos = start_pos
    read_current_pos = 1
    del_len = 0
    for idx, item in enumerate(re.findall(comp, cols[5])):
        span, letter = item[:-1], item[-1]
        if letter in 'MND':
            current_pos += int(span)
        if letter in 'MIS':
            read_current_pos += int(span)
        if letter == 'D':
            del_len += int(span)
        if letter == 'N' and int(span) > 20000 and 10 < start_pos < 150:
            gap_len = int(span)
            read_current_index = read_current_pos-1
            current_index = current_pos-1
            near_query_seq = cols[9][read_current_index-min(30, read_current_index):read_current_index+kmer+30]
            near_genome_region = seq[current_index-gap_len-30:current_index-gap_len]+seq[current_index:current_index+kmer+30]
            if near_query_seq:
                print('\t'.join([cols[0], cols[1], cols[3], str(current_pos), near_query_seq, near_genome_region]))
                new_cols = cols
                new_cols[3] = str(start_pos + gap_len)
                new_cols[5] = ''.join([i for i in re.findall(comp, cols[5]) if 'N' not in i])
                g.write('\t'.join(new_cols[:11]+[str(current_pos)])+'\n')
                h.write(line)
            break
f.close()
g.close()
h.close()