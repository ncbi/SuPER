# -*- coding: utf-8 -*-
"""
Created on Tue May 26 21:13:09 2020

@author: LucasTsubasaYang
"""
import re, sys

left_XA_file = sys.argv[1]
right_XA_file = sys.argv[2]
both_SA_file = sys.argv[3]
filter_sam_file = sys.argv[4]
spliced_reads_file = sys.argv[5]

def update_CIGAR_startpos(str1, str2, pos):
    str1 = str1.replace('H', 'S')
    pos = int(pos)
    comp = re.compile('[0-9]+[SHMNDIPX]')
    str1_list = re.findall(comp, str1)
    str2_list = re.findall(comp, str2)
    length = 0
    str2_s_len = int(str2_list[0][:-1])
    
    new_str2 = ''
    for item in str1_list:
        span, letter = item[:-1], item[-1]
        if letter in 'MIS':
            length += int(span)
        if length > str2_s_len:
            diff = length - str2_s_len
            if letter == 'S':
                # if still got S, stop updating and only keep the right read
                return str2, str(pos)
            new_overlap_str = str(int(span)-diff) + letter
#            print(span, diff, str2_s_len, pos)
            new_str2 += new_overlap_str
            # update startpos
            for i in re.findall(comp, new_str2):
                s, l = i[:-1], i[-1]
                if l in 'MND':
                    pos -= int(s)
            break
        new_str2 += item
    new_str2 += "".join(str2_list[1:])
    
    return new_str2, str(pos)


g = open(filter_sam_file, 'w') # $12 must be the smuggling split site
h = open(spliced_reads_file, 'w') # $4 pos must be original one

# left XA file
f = open(left_XA_file)
for line in f:
    cols = line.strip().split('\t')
    suppl_cigar, suppl_pos = cols[-2], cols[-1] # right cigar
    cigar = cols[5] # left cigar
    cols[3] = suppl_pos # right pos
    pos = cols[3]
    res = update_CIGAR_startpos(cigar, suppl_cigar, pos)
    if res:
        h.write('\t'.join(cols)+'\n')
        new_cigar, new_pos = res
        cols[3] = new_pos
        cols[5] = new_cigar
        g.write('\t'.join(cols[:11]+[pos])+'\n')
f.close()

# right XA file
f = open(right_XA_file)
for line in f:
    cols = line.strip().split('\t')
    suppl_cigar = cols[-2] # left cigar
    cigar = cols[5] # right cigar
    pos = cols[3] # right pos
    res = update_CIGAR_startpos(suppl_cigar, cigar, pos)
    if res:
        h.write('\t'.join(cols)+'\n')
        new_cigar, new_pos = res
        cols[3] = new_pos
        cols[5] = new_cigar
        g.write('\t'.join(cols[:11]+[pos])+'\n')
f.close()

# both SA file
f = open(both_SA_file)
for line in f:
    cols = line.strip().split('\t')
    suppl_cigar, suppl_pos = cols[-2], cols[-1] # right cigar
    cigar = cols[5] # left cigar
    cols[3] = suppl_pos # right pos
    pos = cols[3]
    res = update_CIGAR_startpos(cigar, suppl_cigar, pos)
    if res:
        h.write('\t'.join(cols)+'\n')
        new_cigar, new_pos = res
        cols[3] = new_pos
        cols[5] = new_cigar
        g.write('\t'.join(cols[:11]+[pos])+'\n')
f.close()

g.close()
h.close()