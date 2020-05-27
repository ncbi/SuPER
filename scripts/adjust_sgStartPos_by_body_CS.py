#-*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:40:54 2020

@author: LucasTsubasaYang
"""

from Bio import SeqIO
import sys, distance

filtered_reads_path = sys.argv[1]
leaderCS = sys.argv[2]
ref_path = sys.argv[3]

kmer = len(leaderCS)

for rc in SeqIO.parse(ref_path, 'fasta'):
    seq = str(rc.seq)

def find_putative_CS(start, end, kmer, genome, leader_core_sequence):
    search_seq = genome[start:end]
    min_dist = kmer
    res_seq = ''
    res_index = 0
    find = 0
    for i in range(0, len(search_seq)-kmer+1):
        query_seq = search_seq[i:i+kmer]
        index = start + i + 1
        dist = distance.hamming(leader_core_sequence, query_seq)
        if dist <= min_dist or leaderCS[-4:] == query_seq[-4:]:
            min_dist = dist
            res_seq = query_seq
            res_index = index
    if min_dist <= 2:
        find = 1
    return res_seq, min_dist, res_index, find

# filtered_reads_path has 2 columns: readName, startPos
subgenome_index_dict = {}
f = open(filtered_reads_path)
for line in f:
    cols = line.strip().split()
    read_name, sg_pos = cols
    if sg_pos not in subgenome_index_dict:
        subgenome_index_dict[sg_pos] = [read_name]
    else:
        subgenome_index_dict[sg_pos].append(read_name)

print('\t'.join(['Adjusted_subgenome_position', 'Original_subgenome_position', 'Putative_TRS-B', 'Hamdist_from_leader_CS', 'Supported_read_list']))
for sg_pos in subgenome_index_dict:
    index = int(sg_pos) - 1
    start = index-15
    end = index+kmer+20
    read_list = subgenome_index_dict[sg_pos]
    csb, hamdist, adjust_sg_pos, find = find_putative_CS(start, end, kmer, seq, leaderCS)
    if find: 
        print('\t'.join([str(adjust_sg_pos), str(sg_pos), csb, str(hamdist), ','.join(read_list)]))
    else:
        print('\t'.join([str(sg_pos), str(sg_pos), 'Not found', 'NA', ','.join(read_list)]))

