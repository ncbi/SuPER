# -*- coding: utf-8 -*-
"""
Created on Fri May 15 16:51:52 2020

@author: LucasTsubasaYang
"""
import os, sys, subprocess, time
from Bio import SeqIO

def cmalign(fasfile, cmfile, outfile):
    cmd = 'cmalign '+ cmfile+' '+fasfile+' > '+outfile
    p = subprocess.Popen(cmd, bufsize=-1, shell=True, universal_newlines=True, 
                         stdout=subprocess.PIPE, executable='/bin/bash')
    output = p.communicate()[0]
    return output

def read_cmalign_result(outfile):
    seq_dict = {}
    f = open(outfile)
    for line in f:
        if line.startswith('//'):
            break
        if line.startswith('#'):
            if 'GC RF' in line:
                name, label, seq = line.strip().split()
                name = name + ' ' + label
                if name in seq_dict:
                    seq_dict[name] += seq
                else:
                    seq_dict[name] = seq
            else:
                continue
        elif line.strip()!='':
            name, seq = line.strip().split()
            if name in seq_dict:
                seq_dict[name] += seq
            else:
                seq_dict[name] = seq
    return seq_dict

def get_TRS_L(cmalign_file, ID, grp, refseq):
    grp_TRS_L_dict = {'Alpha': 'AACUAAAC', 'Beta': 'AACG.AAC', 
                      'Gamma': 'aUUAAaa', 'Delta': 'GACACC.A'}
    grp_TRS_L = grp_TRS_L_dict[grp]
    seq_dict = read_cmalign_result(cmalign_file)
    try:
        TRS_L_index = seq_dict['#=GC RF'].index(grp_TRS_L)
    except:
        grp_TRS_L = grp_TRS_L.replace('.','')
        TRS_L_index = seq_dict['#=GC RF'].index(grp_TRS_L)
    new_TRS_L = seq_dict[ID][TRS_L_index:TRS_L_index+len(grp_TRS_L)]
    new_TRS_L = new_TRS_L.replace('U','T').upper()
    try:
        position = str(refseq.index(new_TRS_L)+1)
    except:
        position = 'Not found'
    return grp_TRS_L, new_TRS_L, position

def main(infile, grp, prefix, cm_dir, firstbp=300):
    cmfile = os.path.join(cm_dir, grp+'_5UTR.cm')
    outfile = prefix+'.cmaln.txt'
    
    # extract 5UTR region
    rc = SeqIO.read(infile, 'fasta')
    refseq = str(rc.seq)
    ID = rc.id
    fiveUTR_region = rc[:firstbp]
    fasfile = prefix+'.5UTR.fas'
    SeqIO.write(fiveUTR_region, fasfile, 'fasta')

    # cmalign
    out = cmalign(fasfile, cmfile, outfile)
    print(out)
    time.sleep(1)
    
    # interpret cmalign result
    grp_TRS_L, new_TRS_L, position = get_TRS_L(outfile, ID, grp, refseq)
    return ID, grp, new_TRS_L, grp_TRS_L, position

if __name__ == '__main__':
    # infile: genome fasta file
    # grp: the group of coronavirus could only be 'Alpha', 'Beta', 'Gamma', 'Delta'
    # outfile: file to store TRS-L and its position
    infile = sys.argv[1]
    grp = sys.argv[2]
    prefix = sys.argv[3]
    outfile = sys.argv[4]
    
    work_dir = os.path.dirname(sys.argv[0]).replace('/scripts','')
    cm_dir = work_dir+'/Rfam_CM/'
    ID, grp, csl, grp_TRS_L, position = main(infile, grp, prefix, cm_dir)
    g = open(outfile, 'w')
    g.write("#TRS-L: "+csl+"\n"+"#Position: "+position+"\n")
    g.close()