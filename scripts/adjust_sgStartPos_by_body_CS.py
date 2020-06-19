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
        if dist <= min_dist:
            min_dist = dist
            res_seq = query_seq
            res_index = index
    if min_dist <= 2:
        find = 1
    return res_seq, str(min_dist), str(res_index), find

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

print('\t'.join(['Adjusted_TRS-B_position', 'Original_TRS-B_position', 'Putative_TRS-B', 'Hamdist_from_TRS-L', 'Supported_read_list']))
adjust_sg_pos_dict = {}
for sg_pos in sorted(subgenome_index_dict.keys()):
    index = int(sg_pos) - 1
    start = index-30
    end = index+kmer+30
    read_list = subgenome_index_dict[sg_pos]
    csb, hamdist, adjust_sg_pos, find = find_putative_CS(start, end, kmer, seq, leaderCS)
    sg_pos = str(sg_pos)
    
    if find:
        # adjust_sg_pos_dict is used to remove site redundancy
        # if adjust_sg_pos is the same, just keep the first original sg_pos and combine the read_list
        if adjust_sg_pos not in adjust_sg_pos_dict:
            adjust_sg_pos_dict[adjust_sg_pos] = [adjust_sg_pos, sg_pos, csb, hamdist, ','.join(read_list)]
        else:
            adjust_sg_pos_dict[adjust_sg_pos][4] += ','.join(read_list)
    else:
        # if orf is not found and the new sg_pos is within 10bp of exisiting sg_pos, its read_list will be acquired by the existing one
        for sg_pos_i in adjust_sg_pos_dict:
            if int(sg_pos_i)-10 <= int(sg_pos) <= int(sg_pos_i)+10:
                adjust_sg_pos_dict[sg_pos_i][4] += ','.join(read_list)
                break
        # the new sg_pos is not within 10bp of any old sg_pos
        else:
            adjust_sg_pos_dict[sg_pos] = [sg_pos, sg_pos, 'Not found', 'NA', ','.join(read_list)]
        
for adjust_sg_pos in adjust_sg_pos_dict:
    res_list = adjust_sg_pos_dict[adjust_sg_pos]
    print('\t'.join(res_list))
    