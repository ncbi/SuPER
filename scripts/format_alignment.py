#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os, sys
from Bio import pairwise2
from Bio import SeqIO
from Bio.Alphabet import IUPAC

TRS_L_file = sys.argv[1]
sgmRNA_doc = sys.argv[2]
ref_file = sys.argv[3]

def locate_pos(align1, align2, score, begin, end): 
    start1 = len(align1[:begin]) - align1[:begin].count("-")
    start2 = len(align2[:begin]) - align2[:begin].count("-")
    return[start1,start1+end-begin,start2,start2+end-begin]

def format_align(aln1,aln2):
    as1,ae1,as2,ae2 = locate_pos(*aln1)
    bs1,be1,bs2,be2 = locate_pos(*aln2)

    lead = aln1[0].strip("-").upper()
    sgmRNA1 = aln1[1].strip("-").upper()
    sgmRNA2 = aln2[0].strip("-").upper()
    body  = aln2[1].strip("-").upper()

    start_m = 12 # max(len(str(as1)), len(str(bs2)), 6) + 1

    if as2 >= bs1:
        l_blank=as2-bs1
        b_blank=0
    else:
        l_blank=0
        b_blank=bs1-as2

    s1_line = ["{:>{width}}".format("TRS-L "+str(as1+1)+" ", width=start_m+l_blank)]  # seq1 line 
    m_line1 = [" " * (start_m+l_blank)]  # match line 
    s2_line = ["{:>{width}}".format("sgmRNA ", width=start_m)]  # seq2 line 
    s1_line += lead[as1:ae1]
    s2_line += sgmRNA1[min(as2,bs1):max(ae2,be1)]

    for n, (a, b) in enumerate(zip(lead[as1:ae1], 
                                   sgmRNA1[as2:ae2])): 
        # Since list elements can be of different length, we center them, 
        # using the maximum length of the two compared elements as width 
        m_len = max(len(a), len(b)) 
        if a == b: 
            m_line1.append("{:^{width}}".format("|", width=m_len))  # match 
        elif a.strip() == "-" or b.strip() == "-": 
            m_line1.append("{:^{width}}".format(" ", width=m_len))  # gap 
        else: 
            m_line1.append("{:^{width}}".format(".", width=m_len))  # mismatch 
            
    first="\n".join(["".join(s1_line), "".join(m_line1), "".join(s2_line)]) +"\n"

    s3_line = ["{:>{width}}".format("TRS-B "+str(bs2+1)+" ", width=start_m+b_blank)]  # seq2 line 
    m_line2 = [" " * (start_m+b_blank)]
#    s4_line = ["{:>{width}}".format("sgmRNA ", width=start_m)]  # seq1 line
    
    s3_line += body[bs2:be2]
#    s4_line += sgmRNA2[min(as2,bs1):max(ae2, be1)]

    for n, (a, b) in enumerate(zip(sgmRNA2[bs1:be1], 
                                   body[bs2:be2])):
        # Since list elements can be of different length, we center them, 
        # using the maximum length of the two compared elements as width 
        m_len = max(len(a), len(b)) 
        if a == b: 
            m_line2.append("{:^{width}}".format("|", width=m_len))  # match 
        elif a.strip() == "-" or b.strip() == "-": 
            m_line2.append("{:^{width}}".format(" ", width=m_len))  # gap 
        else: 
            m_line2.append("{:^{width}}".format(".", width=m_len))  # mismatch 

    return first+"\n".join(["".join(m_line2), "".join(s3_line)])

# read leading sequence
seq1 = SeqIO.read(TRS_L_file, "fasta", alphabet=IUPAC.ambiguous_dna)
# read the whole genome
seq3= SeqIO.read(ref_file, "fasta", alphabet=IUPAC.ambiguous_dna)
# read sgmRNA (more region match the TRS-B than TRS-L)
print("="*36+"Subgenomic mRNA Alignment"+"="*36)
for fn in sorted(os.listdir(sgmRNA_doc)):
    sgmRNA_file = os.path.join(sgmRNA_doc, fn)
    seq2 = SeqIO.read(sgmRNA_file, "fasta", alphabet=IUPAC.ambiguous_dna)
    
    aln1=pairwise2.align.localms(seq1.seq,seq2.seq,3,-2,-50,-10)
    aln2=pairwise2.align.localms(seq2.seq,seq3.seq,3,-2,-50,-10)
    
    pos = os.path.basename(fn).split('__')[-1].replace('.fasta','')
    print("Position "+pos+" alignment:")
    print(format_align(aln1=aln1[0],aln2=aln2[0]))
    print()
    print("="*100)
