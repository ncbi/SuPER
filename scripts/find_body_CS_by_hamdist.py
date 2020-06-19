# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 14:19:21 2020

@author: LucasTsubasaYang
"""

from Bio import SeqIO
import sys, distance

leaderCS = sys.argv[1]
fas_path = sys.argv[2]
gff_path = sys.argv[3]
kmer = len(leaderCS)

for rc in SeqIO.parse(fas_path, 'fasta'):
    seq = str(rc.seq)

if gff_path == 'None':
    # without ORF information
    print('\t'.join(['Adjusted_TRS-B_position','Original_TRS-B_position','Putative_TRS-B','Hamdist_from_TRS-L']))
    for i in range(0, len(seq)-kmer+1):
        if i > 20000 and i < len(seq)-90:
            subseq = seq[i:i+kmer]
            dist = distance.hamming(leaderCS, subseq)
            if dist <= 2:
                print('\t'.join([str(i+1),str(i+1),subseq,str(dist)]))
else:
    # with ORF information
    # read gff file
    f = open(gff_path)
    gene_range = {}
    for line in f:
        if line.strip()!='' and not line.startswith('#'):
            cols = line.strip().split('\t')
            if cols[2] == 'CDS':
                start, end = cols[3], cols[4]
                if 'product=' in cols[8]:
                    gene_name =[kv.split("=") for kv in cols[8].split(';') if kv.startswith('product=')][0][-1]
                elif 'gene=' in cols[8]:
                    gene_name = [kv.split("=") for kv in cols[8].split(';') if kv.startswith('gene=')][0][-1]
                else:
                    gene_name = 'unknown_gene'
                gene_range[gene_name+":"+start+"-"+end] = (start, end)
    f.close()

    print('\t'.join(['Adjusted_TRS-B_position','Original_TRS-B_position','Putative_TRS-B','Hamdist_from_TRS-L','Recommend_label*','Distance_from_ORF','Associated_ORF:start-end']))
    sg_pos_dict = {}
    for i in range(0, len(seq)-kmer+1):
        if i > 20000 and i < len(seq)-90:
            find_ORF = 0
            subseq = seq[i:i+kmer]
            position = i+1
            dist = distance.hamming(leaderCS, subseq)
            if dist <= 2:
                # with ORF support the max distance can be 2
                for gene_name in gene_range:
                    start, end = gene_range[gene_name]
                    if 0 < int(start)-position < 170:
                        find_ORF = 1
                        orf_startpos = gene_name.split(":")[-1].split("-")[0]
                        dist_orf = int(orf_startpos)-position
                        if dist > 1:
                            recomd_lv = 'Not recommended'
                        else:
                            recomd_lv = '-'
                        sg_pos_dict[str(position)] = [str(position),str(position),subseq,str(dist), recomd_lv, str(dist_orf), gene_name]
                        break
                else:
                    # without ORF support the max distance can only be 1
                    if dist <= 1:
                        # because dist <= 1 the recomd_lv is always '-'
                        sg_pos_dict[str(position)] = [str(position),str(position),subseq,str(dist), "-", "NA", "Not found"]

    # remove sg_pos redundancy, if the orf is the same then we get the min dist pos
    orf_dict = {}
    for sg_pos in sg_pos_dict:
        sg_pos, sg_pos, trsb, dist, dist_orf, recomd_lv, orf = sg_pos_dict[sg_pos]
        if orf == "Not found":
#            orf_dict[orf+"_"+sg_pos] = [sg_pos, sg_pos, trsb, dist, orf] # no output for those without corresponding orf
            continue
        if orf not in orf_dict:
            orf_dict[orf] = [sg_pos, sg_pos, trsb, dist, dist_orf, recomd_lv, orf]
        else:
            if int(dist) < int(orf_dict[orf][3]):# update value to get smaller distance TRS-B
                orf_dict[orf] = [sg_pos, sg_pos, trsb, dist, dist_orf, recomd_lv, orf]
            elif int(dist) == int(orf_dict[orf][3]):
                if int(sg_pos) > int(orf_dict[orf][1]):# when with same distance, update value to make TRS-B closer to ORF
                    orf_dict[orf] = [sg_pos, sg_pos, trsb, dist, dist_orf, recomd_lv, orf]
    
    for orf in orf_dict:
        res_list = orf_dict[orf]
        print('\t'.join(res_list))
    print("*Notice: The TRS-B with hamming distance from TRS-L greater than 1 would be labeled as 'Not recommended'.")
