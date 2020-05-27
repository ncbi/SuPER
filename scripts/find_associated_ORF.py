# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:45:19 2020

@author: LucasTsubasaYang
"""
import sys

# file including subgenome positions found by RNAseq or Nanopore
subgenome_pos_file = sys.argv[1]
# genome annotation gff file
gff_path = sys.argv[2]

f = open(gff_path)
gene_range = {}
for line in f:
    if line.strip()!='' and not line.startswith('#'):
        cols = line.strip().split('\t')
        if cols[2] == 'CDS':
            start, end = cols[3], cols[4]
            try:
                gene_name = [kv.split("=") for kv in cols[8].split(';') if kv.startswith('gene=')][0][-1]
            except:
                gene_name =[kv.split("=") for kv in cols[8].split(';') if kv.startswith('product=')][0][-1]
            gene_range[gene_name+':'+start+'-'+end] = (start, end)
f.close()

f = open(subgenome_pos_file)
line = f.readline()
cols = line.strip().split('\t')
print('\t'.join(cols[:-1]+['Associated_ORF:start-end', cols[-1]]))
for line in f:
    cols = line.strip().split('\t')
    adjust_sg_pos = cols[0]
    find_ORF = 0
    for gene_name in gene_range:
        start, end = gene_range[gene_name]
        if 0 < int(start)-int(adjust_sg_pos) < 100:
            print('\t'.join(cols[:-1]+[gene_name, cols[-1]]))
            find_ORF = 1
            break
    if not find_ORF:
        print('\t'.join(cols[:-1]+['Not found', cols[-1]]))
f.close()
