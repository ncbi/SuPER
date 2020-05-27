# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 16:08:25 2020

@author: LucasTsubasaYang
"""
import sys, re, math

def update_coverage(start, end, pos_cov_dict):
    for i in range(start, end+1):
        pos_cov_dict[i] += 1
    return pos_cov_dict

filtered_sam_path = sys.argv[1]
cov_outfile = sys.argv[2]
ref_path = sys.argv[3]
read_sum_num = int(sys.argv[4])

seq = ''
f = open(ref_path)
for line in f:
    if not line.startswith('>'):
        seq += line.strip()
genome_len = len(seq)

# initialize all pos as keys in pos_cov_dict
pos_cov_dict = {}
for i in range(1, genome_len+1):
    pos_cov_dict[i] = 0

# only focus on MND
comp = re.compile('[0-9]+[MND]')
split_reads_supported_sites = {}
f = open(filtered_sam_path)
for line in f:
    if line.startswith('@'):
        continue
    cols = line.strip().split('\t')
    start_pos = int(cols[3])
    end_pos = start_pos-1
    for item in re.findall(comp, cols[5]):
        length, letter = item[:-1], item[-1]
        if letter != 'N':
            start_pos = end_pos+1
            end_pos += int(length)
            pos_cov_dict = update_coverage(start_pos, end_pos, pos_cov_dict)
        else:
            start_pos = end_pos+1
            end_pos += int(length)
f.close()

pos_list = []
pre_log_cov = 0
cov_cutoff = 1
g = open(cov_outfile,'w')
for pos, cov in sorted(pos_cov_dict.items(),key=lambda item:item[0], reverse=False):
    cov = pos_cov_dict[pos]
    try:
        log_cov = math.log(cov, 10)
    except:
        log_cov = 0
    if pos > 20000 and cov/read_sum_num > 0:
        cov_cutoff = cov/read_sum_num
    if pos > 20000 and log_cov-pre_log_cov > 0.03 and cov/read_sum_num >= cov_cutoff:
        pos_list.append(pos)
        print('\t'.join([str(pos), str(cov), str(log_cov)]))
    g.write('\t'.join([str(pos), str(cov), str(log_cov)])+'\n')
    pre_log_cov = log_cov
g.close()
