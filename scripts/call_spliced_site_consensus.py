# -*- coding: utf-8 -*-
"""
Created on Wed May 13 13:43:43 2020

@author: LucasTsubasaYang
"""
import sys, os, shutil
import subprocess

def run_cmd(cmd):
    p = subprocess.Popen(cmd,bufsize=-1, shell=True, universal_newlines=True, stdout=subprocess.PIPE,executable='/bin/bash')
    stdout, stderr = p.communicate()
    code = p.returncode
    if code != 0:
        raise SystemExit('Error: {0}'.format(stderr))

if __name__ == "__main__":
    # parameters
    sam_file = sys.argv[1]
    header = sys.argv[2]
    pos_list = sys.argv[3].split(',')
    prefix = sys.argv[4]
    type_prefix = sys.argv[5]
    py_inter = sys.argv[6]
    
    script_dir = os.path.dirname(sys.argv[0])
    tmpdir = os.path.dirname(prefix)
    filtered_sam_dir = tmpdir+'/filtered_sam'
    consensus_dir = tmpdir+'/consensus'
    if os.path.exists(filtered_sam_dir) and os.listdir(filtered_sam_dir)!=[]:
        shutil.rmtree(filtered_sam_dir)
    if os.path.exists(consensus_dir) and os.listdir(consensus_dir)!=[]:
        shutil.rmtree(consensus_dir)

    # step.1: get non-redundant sites
    pos_dict = {}
    for pos in pos_list:
        pos = int(pos)
        for ex_pos in pos_dict:
            if abs(pos-ex_pos) < 10:
                pos_dict[ex_pos].append(pos)
                break
        else:
            pos_dict[pos] = [pos]

    # step.2: get valid reads according to per position
    # Notice: the sam header must write so the sam2consensus can read.
    for key in pos_dict:
        print('Dealing with pos', key)
        # $12 is a smuggling information to record the split position of the read
        if type_prefix == 'nanopore':
            # nanopore data the $12 must be equal to the position
            cmd = '''
            mkdir -p {tmpdir}/filtered_sam
            cat {header} > {tmpdir}/filtered_sam/{key}.filtered.sam
            awk 'BEGIN{{FS="\\t";OFS="\\t"}}$12=={key}{{print $0}}' {sam_file} >> {tmpdir}/filtered_sam/{key}.filtered.sam
            '''.format(key=key, sam_file=sam_file, header=header, prefix=prefix, tmpdir=tmpdir, type_prefix=type_prefix)
            print(cmd)
            run_cmd(cmd)
        else:
            # other data type can be in a range
            cmd = '''
            mkdir -p {tmpdir}/filtered_sam
            cat {header} > {tmpdir}/filtered_sam/{key}.filtered.sam
            awk 'BEGIN{{FS="\\t";OFS="\\t"}}$12<=({key}+10)&&$12>=({key}-10){{print $0}}' {sam_file} >> {tmpdir}/filtered_sam/{key}.filtered.sam
            '''.format(key=key, sam_file=sam_file, header=header, prefix=prefix, tmpdir=tmpdir, type_prefix=type_prefix)
            print(cmd)
            run_cmd(cmd)

     # step.3: find consensus
        if type_prefix == 'nanopore':
            cutoff = 0.50 # due to the high error rate in nanopore data so we set it low
            cmd = '''
            mkdir -p {tmpdir}/consensus
            {py_inter} {script_dir}/sam2consensus.py -i {tmpdir}/filtered_sam/{key}.filtered.sam -c {cutoff} -p {key} -o {tmpdir}/consensus
            '''.format(py_inter=py_inter, script_dir=script_dir, key=key, prefix=prefix, type_prefix=type_prefix, tmpdir=tmpdir, cutoff=cutoff)
            print(cmd)
            run_cmd(cmd)
            
            # extract the split sites nearby region from long ONT data
            path = [os.path.join(tmpdir+'/consensus',fn) for fn in os.listdir(tmpdir+'/consensus') if '__'+str(key) in fn][0]
            f = open(path)
            lines = f.readlines()
            Id_line = lines[0]
            seq = lines[1].strip()
            subseq = seq[key-1-40:key-1+60]
            f.close()
            g = open(path, 'w')
            g.write(Id_line)
            g.write(subseq.upper()+'\n')
            g.close()
            
        else:
            cutoff = 0.75
            cmd = '''
            mkdir -p {tmpdir}/consensus
            {py_inter} {script_dir}/sam2consensus.py -i {tmpdir}/filtered_sam/{key}.filtered.sam -c {cutoff} -f '' -p {key} -o {tmpdir}/consensus
            '''.format(py_inter=py_inter, script_dir=script_dir, key=key, prefix=prefix, type_prefix=type_prefix, tmpdir=tmpdir, cutoff=cutoff)
            print(cmd)
            run_cmd(cmd)
