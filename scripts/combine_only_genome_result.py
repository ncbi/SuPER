# -*- coding: utf-8 -*-
"""
Created on Fri Jun 12 09:45:13 2020

@author: LucasTsubasaYang
"""

"""This program is used to combine positions found by RNA-seq data and sequence similarity"""
import sys

rnaseq_file = sys.argv[1]
seqsim_file = sys.argv[2]


header = ['Adjusted_TRS-B_position','Original_TRS-B_position','Putative_TRS-B','Hamdist_from_TRS-L', 'Recommend_label*', \
          'Distance_from_ORF', 'Associated_ORF:start-end','Supported_by_RNAseq','Supported_read_list']

pos_list = []
seqsim_dict = {}
f1 = open(seqsim_file)
f1.readline() # pass header
for line in f1:
    if line.startswith("*"):
        continue
    cols = line.strip().split('\t')
    pos, pos_prime, trsb, ham_dist, rec_lv, dist_orf, orf  = cols
    if orf != 'Not found':
        pos_list.append(pos)
        seqsim_dict[pos] = [pos, pos_prime, trsb, ham_dist, rec_lv, dist_orf, orf, 'No', 'NA']
f1.close()

rnaseq_dict = {}
f2 = open(rnaseq_file) # pass header
f2.readline()
for line in f2:
    cols = line.strip().split('\t')
    adjust_pos, original_pos, trsb, ham_dist, dist_orf, orf, read_lst = cols
    if ham_dist == 'NA' or int(ham_dist) > 1:
        recomd_lv = 'Not recommended'
    else:
        recomd_lv = '-'
    if adjust_pos not in pos_list:
        pos_list.append(adjust_pos)
    rnaseq_dict[adjust_pos] = [adjust_pos, original_pos, trsb, ham_dist, recomd_lv, dist_orf, orf, 'Yes', read_lst]
f2.close()

print('\t'.join(header))
for pos in sorted(pos_list):
    if pos in rnaseq_dict:
        print('\t'.join(rnaseq_dict[pos]))
    else:
        print('\t'.join(seqsim_dict[pos]))
print("*Notice: The TRS-B with hamming distance from TRS-L greater than 1 would be labeled as 'Not recommended'.")
