# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:39:20 2020

@author: LucasTsubasaYang
"""
import os, sys, subprocess, tempfile, time
from Bio import SeqIO
from distutils.spawn import find_executable
from progress.spinner import Spinner

def run_cmd(cmd, wait=False):
    p = subprocess.Popen(cmd,bufsize=-1, shell=True, universal_newlines=True, stdout=subprocess.PIPE, 
                         stderr=subprocess.PIPE, executable='/bin/bash')
    if wait:
        spin = Spinner("Processing...")
        while p.poll() is None:
            time.sleep(0.2)
            spin.next()
        print()
    stdout, stderr = p.communicate()
    code = p.returncode
    if code != 0:
        raise SystemExit('Error: {0}'.format(stderr))

def verbose_print(msg, cmd, verbose):
    if verbose:
        print(msg+' '+cmd)
    else:
        print(msg+' running...')

def is_tool(name):
    """Check whether program 'name' is on PATH."""
    return find_executable(name) is not None

def pipeline(args):
    inputfile=args.input
    ref=args.reference
    grp=args.group
    datatype=args.datatype
    gff=args.gff
    csl=args.csl
    program=args.program
    cutoff=args.cutoff
    verbose=args.verbose
    core=args.core
    tmpdir=args.tmpdir
    outputfile=args.output
    
    work_dir = os.path.dirname(sys.argv[0])
    if work_dir == '':
        script_dir = "./scripts"
    else:
        script_dir = work_dir+"/scripts"
        
    f = next(tempfile._get_candidate_names())
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    prefix=tmpdir+"/"+f
    
    print("="*80)
    
    # find python interpreter
    if find_executable('pypy3') is not None:
        py_inter = 'pypy3'
    else:
        if find_executable('python') is not None:
            py_inter = 'python'
        else:
            py_inter = 'python3'
    print("Using Python Interpreter:", py_inter)
    print("="*80)

    if grp not in ['Alpha', 'Beta', 'Gamma', 'Delta']:
        print("[Error]: -g/--group must be Alpha, Beta, Gamma or Delta. Please try again.")
        return
    
    # step.1: infer TRS-L
    rc = SeqIO.read(ref, 'fasta')
    refseq = str(rc.seq)
    if not csl:
        # TRS-L isn't provided by user, then auto detect
        cmd = '''{py_inter} {script_dir}/infer_TRS_L.py {ref_file} {group} {prefix} {prefix}.TRS_L.txt
        '''.format(py_inter=py_inter, script_dir=script_dir, ref_file=ref, group=grp, prefix=prefix)
        if verbose:
            msg = cmd
        else:
            msg = 'Running...'
        print("[Infer TRS-L]: ",msg)
        run_cmd(cmd)
        lines = open(prefix+'.TRS_L.txt').readlines()
        csl = lines[0].strip().split(": ")[-1]
        position = lines[1].strip().split(": ")[-1]
        if '-' in csl:
            print("[Error]: Detected putative TRS-L isn't correct! Please offer one by argument -l/--csl.")
            return
        else:
            print("[Infer TRS-L]: Detect",csl,"as TRS-L. If you suggest another one, please offer it by argument -l/--csl.")
    else:
        try:
            position = str(refseq.index(csl)+1)
        except:
            position = 'Not found'
        g = open('{prefix}.TRS_L.txt'.format(prefix=prefix), 'w')
        g.write("#TRS-L: "+csl+"\n"+"#Position: "+position+"\n")
        g.close()
        print("[Infer TRS-L]:", csl, "will be used as TRS-L.")
    kmer = len(csl)
    print("[Infer TRS-L]: done!")
    print("="*80)

    # step.2.1: if SAM file is not given, search TRS-B through the genome
    if not inputfile:
        if gff:
            gff_path = gff
        else:
            gff_path = 'None'
        cmd = '''
        {py_inter} {script_dir}/find_body_CS_by_hamdist.py {CS_L} {ref_file} {gff_path} > {prefix}.genome_subgenome_pos.tab
        cat {prefix}.TRS_L.txt {prefix}.genome_subgenome_pos.tab > {output}
        '''.format(py_inter=py_inter, script_dir=script_dir, CS_L=csl, ref_file=ref, gff_path=gff_path, prefix=prefix, output=outputfile)
        if verbose:
            msg = cmd
        else:
            msg = 'Running...'
        print("[(Only genome) Find TRS-B positions]:", msg)
        run_cmd(cmd)
        print("[(Only genome) Find TRS-B positions]: done!")
        
    else:
    # step.2.1: firstly search TRS-B through the genome
        if gff:
            gff_path = gff
            cmd = '''
            {py_inter} {script_dir}/find_body_CS_by_hamdist.py {CS_L} {ref_file} {gff_path} > {prefix}.genome_subgenome_pos.tab
            '''.format(py_inter=py_inter, script_dir=script_dir, CS_L=csl, ref_file=ref, gff_path=gff_path, prefix=prefix, output=outputfile)
            if verbose:
                msg = cmd
            else:
                msg = 'Running...'
            print("[Find TRS-B positions]:", msg)
            run_cmd(cmd)
            print("[Find TRS-B positions]: done!")
            
    # step.2.2: infer align program, different program will give diff modes of reads
    #         for example: bwa tends to give partially mapped reads while hisat2 and minimap2 
    #         tend to add a large span of gaps into the reads (i.e. junctions)
        if not program:
            # if not given then auto detect
            pg_line = ''
            f = open(inputfile)
            for line in f:
                if line.startswith("@PG"):
                    pg_line = line.strip()
                    break
            f.close()
            program = pg_line.split('\t')[1].replace('ID:','')
            if len(pg_line) == 0:
                print("[Error]: Sam file without header! or alignment program is not found in sam file! Please offer it by argument -g or --program.")
                return
        print("[Infer alignment program]: Use", program,"as alignment program.")
        print("[Infer alignment program]: done!")
        print("="*80)
        
        # step.3: deal with alignment data to get header
        sam_processed_file = inputfile
        cmd = '''head -n10 {sam_file} | grep -e "^@" > {prefix}.sam.header'''.format(sam_file=inputfile, prefix=prefix, core=core)
        if verbose:
            msg = cmd
        else:
            msg = 'Running...'
        print("[Preprocess alignment file]:",msg)
        run_cmd(cmd)
        print("[Preprocess alignment file]: done!")
        print("="*80)
        
        if datatype:
            type_prefix = 'rnaseq'
            if program in ['bwa', 'bowtie2', 'bowtie']:
            # step.3.1.2: deal with illumina RNA-seq data by bwa/bowtie
                # TRS-B reads
                # get spliced sites reads distribution
                ref_len = len(str(SeqIO.read(ref, 'fasta').seq))
                cmd = '''awk '$3>0 && ($4>20000 && $4<({ref_len}-90)) && ($6~/^([0-9]+)S/ || $6~/^([0-9]+)H/){{print $0}}' {sam_file} | grep -E "XA:Z|SA:Z" > {prefix}.right.sam
                awk '$3>0 && ($4<150 && $4>10) && ($6~/S$/ || $6~/H$/){{print $0}}' {sam_file} | grep -E "XA:Z|SA:Z" > {prefix}.left.sam
                
                # get {prefix}.right spliced read records/sites
                grep "XA:Z" {prefix}.right.sam > {prefix}.right.tmp.sam
                # $3,$2: CIGAR, strand+pos
                grep "XA:Z" {prefix}.right.sam | awk '{{print $NF}}' | awk -F"," 'BEGIN{{OFS="\\t"}}{{print $3,$2}}' > {prefix}.right.tab
                paste {prefix}.right.tmp.sam {prefix}.right.tab | awk 'BEGIN{{OFS="\\t"}}(and(16,$2) && substr($NF,1,1)=="-") || (!and(16,$2) && substr($NF,1,1)=="+"){{$NF=substr($NF,2,length($NF)); print $0}}'| awk '($NF>20000 && $NF<({ref_len}-90)){{print $0}}' > {prefix}.right.XA.sam
                
                # get {prefix}.left spliced read records/sites
                grep "XA:Z" {prefix}.left.sam > {prefix}.left.tmp.sam
                # $3,$2: CIGAR, strand+pos
                grep "XA:Z" {prefix}.left.sam | awk '{{print $NF}}' | awk -F"," 'BEGIN{{OFS="\\t"}}{{print $3,$2}}' > {prefix}.left.tab
                paste {prefix}.left.tmp.sam {prefix}.left.tab | awk 'BEGIN{{OFS="\\t"}}(and(16,$2) && substr($NF,1,1)=="-") || (!and(16,$2) && substr($NF,1,1)=="+"){{$NF=substr($NF,2,length($NF)); print $0}}' | awk '($NF>20000 && $NF<({ref_len}-90)){{print $0}}' > {prefix}.left.XA.sam
                
                # get {prefix}.left & {prefix}.right paired read records/sites
                grep "SA:Z" {prefix}.left.sam > {prefix}.left.1.tmp.sam
                grep "SA:Z" {prefix}.right.sam > {prefix}.right.1.tmp.sam
                # $4,$3,$2: CIGAR, strand, pos
                grep "SA:Z" {prefix}.left.sam | awk '{{print $NF}}' | awk -F"," 'BEGIN{{OFS="\\t"}}{{print $4,$3""$2}}' > {prefix}.left.1.tab
                paste {prefix}.left.1.tmp.sam {prefix}.left.1.tab | awk 'BEGIN{{OFS="\\t"}}(and(16,$2) && substr($NF,1,1)=="-") || (!and(16,$2) && substr($NF,1,1)=="+"){{$NF=substr($NF,2,length($NF)); print $0}}' | awk '($NF>20000 && $NF<({ref_len}-90)){{print $0}}' > {prefix}.both.SA.tmp.sam
                # store the longer sequence
                awk 'BEGIN{{OFS="\\t"}}NR==FNR{{a[$1""$4]=$10;next}}($1""$NF in a){{if(length($10)>=length(a[$1""$NF]))\
                {{print $0}}else{{$10=a[$1""$NF]; print $0}}}}' {prefix}.right.1.tmp.sam {prefix}.both.SA.tmp.sam > {prefix}.both.SA.sam
                
                # update CIGAR and startpos to get  *.filtered.sam
                {py_inter} {script_dir}/update_CIGAR_startpos.py {prefix}.left.XA.sam {prefix}.right.XA.sam {prefix}.both.SA.sam {prefix}.filtered.sam {prefix}.rnaseq_spliced_reads.sam
                '''.format(py_inter=py_inter, script_dir=script_dir, sam_file=sam_processed_file, ref_len=ref_len, prefix=prefix)
                if verbose:
                    msg = cmd
                else:
                    msg = 'Running...'
                print("[(RNA-Seq) Find spliced reads]:", msg)
                run_cmd(cmd, wait=True)
                
                # find TRS-B && TRS-L reads and change the CIGAR and startpos and stach into {prefix}.filtered.sam
                cmd = '''#read_sum=$(wc -l {prefix}.rnaseq_spliced_reads.sam | cut -f1);
                cut -f4 {prefix}.rnaseq_spliced_reads.sam | sort | uniq -c > {prefix}.rnaseq_spliced_reads_distr.tab
                cat {prefix}.rnaseq_spliced_reads_distr.tab | awk '$1>={cutoff}{{print $2}}' > {prefix}.rnaseq_high_cov_reads.list
                
                # recruit other reads to help *.filtered.sam better call consensus
                awk 'BEGIN{{OFS="\\t"}} NR==FNR{{a[$2];next}}NR>FNR{{if(($4 in a) && (match($6,/^([0-9]+)M/,m)||$6 ~/^([0-9]+)H/)){{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$4}}}}' {prefix}.rnaseq_spliced_reads_distr.tab {sam_file} >> {prefix}.filtered.sam
                '''.format(sam_file=sam_processed_file, cutoff=cutoff, prefix=prefix)
                if verbose:
                    msg = cmd
                else:
                    msg = 'Running...'
                print("[(RNA-Seq) Find spliced sites]:", msg)
                run_cmd(cmd, wait=True)
                
                cmd = '''awk 'NR==FNR{{a[$1]=$0;next}}$4 in a{{print $1"."$2"\\t"$4}}' {prefix}.rnaseq_high_cov_reads.list {prefix}.rnaseq_spliced_reads.sam > {prefix}.rnaseq_spliced_reads.filtered.tab'''.format(prefix=prefix)
                if verbose:
                    msg = cmd
                else:
                    msg = 'Running...'
                print("[(RNA-Seq) Filter reads]:", msg)
                run_cmd(cmd)
                
            else:
            # step.3.1.1: deal with illumina RNA-seq data by hisat2 and others
                # find TRS-B reads and change the CIGAR and startpos and stach into {prefix}.filtered.sam
                cmd = '''{py_inter} {script_dir}/find_RNASeq_wide_spliced_reads.py {sam_file} {ref_file} {kmer} {prefix}.filtered.sam > {prefix}.rnaseq_spliced_reads.tab
                '''.format(py_inter=py_inter, script_dir=script_dir, sam_file=sam_processed_file, ref_file=ref, kmer=kmer, prefix=prefix)
                if verbose:
                    msg = cmd
                else:
                    msg = 'Running...'
                print("[(RNA-Seq) Find spliced reads]:", msg)
                run_cmd(cmd, wait=True)
                
                cmd = '''#read_sum=$(wc -l {prefix}.rnaseq_spliced_reads.tab | cut -f1);
                cut -f4 {prefix}.rnaseq_spliced_reads.tab | sort | uniq -c > {prefix}.rnaseq_spliced_reads_distr.tab
                cat {prefix}.rnaseq_spliced_reads_distr.tab | sort -nrk1,1 | awk '$1>={cutoff}{{print $2}}' > {prefix}.rnaseq_high_cov_reads.list
                '''.format(prefix=prefix, cutoff=cutoff)
                if verbose:
                    msg = cmd
                else:
                    msg = 'Running...'
                print("[(RNA-Seq) Find spliced sites]:", msg)
                run_cmd(cmd)
                
                cmd = '''awk 'NR==FNR{{a[$1]=$0;next}}$4 in a{{print $1"."$2"\\t"$4}}' {prefix}.rnaseq_high_cov_reads.list {prefix}.rnaseq_spliced_reads.tab > {prefix}.rnaseq_spliced_reads.filtered.tab'''.format(prefix=prefix)
                if verbose:
                    msg = cmd
                else:
                    msg = 'Running...'
                print("[(RNA-Seq) Filter reads]:", msg)
                run_cmd(cmd)
            
            print("[(RNA-Seq) Find TRS-B positions]: done!")
            print("="*80)
        
        else:
            # step.3.2: deal with oxford Nanopore data
            type_prefix = 'nanopore'
            cmd = '''{py_inter} {script_dir}/find_Nanopore_wide_spliced_reads.py {sam_file} {ref_file} {kmer} {prefix}.filtered.sam {prefix}.nanopore_spliced_reads.sam > {prefix}.nanopore_spliced_reads.tab
            '''.format(py_inter=py_inter, script_dir=script_dir, sam_file=sam_processed_file, ref_file=ref, kmer=kmer, prefix=prefix)
            if verbose:
                msg = cmd
            else:
                msg = 'Running...'
            print("[(Nanopore) Find wide spliced reads]:", msg)
            run_cmd(cmd)
            
            cmd = '''#read_sum=$(wc -l {prefix}.nanopore_spliced_reads.tab | cut -f1);
            cut -f4 {prefix}.nanopore_spliced_reads.tab | sort | uniq -c > {prefix}.nanopore_spliced_reads_distr.tab
            cat {prefix}.nanopore_spliced_reads_distr.tab | sort -nrk1,1 | awk '$1>={cutoff}{{print $2}}' > {prefix}.nanopore_high_cov_reads.list
            read_sum=$(wc -l {prefix}.filtered.sam | cut -f1);
            {py_inter} {script_dir}/find_Nanopore_coverage_jump_point.py {prefix}.nanopore_spliced_reads.sam {prefix}.nanopore_cov_distr.tab {ref_file} $read_sum > {prefix}.nanopore_cov_jump_point.list
            awk 'NR==FNR{{a[$1]=$0;next}}$1 in a{{print $0}}' {prefix}.nanopore_high_cov_reads.list {prefix}.nanopore_cov_jump_point.list > {prefix}.nanopore_high_cov_jump_point.list
            '''.format(py_inter=py_inter, script_dir=script_dir, cutoff=cutoff, ref_file=ref, prefix=prefix)
            if verbose:
                msg = cmd
            else:
                msg = 'Running...'
            print("[(Nanopore) Find coverage jump points]:", msg)
            run_cmd(cmd)
            
            cmd = '''awk 'NR==FNR{{a[$1]=$0;next}}$4 in a{{print $1"."$2"\\t"$4}}' {prefix}.nanopore_high_cov_jump_point.list {prefix}.nanopore_spliced_reads.tab > {prefix}.nanopore_spliced_reads.filtered.tab
            '''.format(py_inter=py_inter, prefix=prefix)
            if verbose:
                msg = cmd
            else:
                msg = 'Running...'
            print("[(Nanopore) Filter reads]:", msg)
            run_cmd(cmd)
            
            print("[(Nanopore) Find TRS-B positions]: done!")
            print("="*80)

        # step.4: find CS-B and use them to adjust the subgenome pos
        cmd ='''{py_inter} {script_dir}/adjust_sgStartPos_by_body_CS.py {prefix}.{type_prefix}_spliced_reads.filtered.tab {CS_L} {ref_file} | sort -k1,1n > {prefix}.{type_prefix}_adjusted_subgenome_pos.tab
        '''.format(py_inter=py_inter, script_dir=script_dir, type_prefix=type_prefix, prefix=prefix, CS_L=csl, ref_file=ref)
        if verbose:
            msg = cmd
        else:
            msg = 'Running...'
        print("[Adjust TRS-B position by body core sequence]:", msg)
        run_cmd(cmd)
        print("[Adjust TRS-B position by body core sequence]: done!")
        print("="*80)

        # step.5: check the subgenomes' associated ORF if given gff file
        if gff:
            cmd = '''{py_inter} {script_dir}/find_associated_ORF.py {prefix}.{type_prefix}_adjusted_subgenome_pos.tab {gff_file} > {prefix}.{type_prefix}_adjusted_subgenome_pos_with_ORF.tab
            '''.format(py_inter=py_inter, script_dir=script_dir, gff_file=gff, prefix=prefix, type_prefix=type_prefix)
            if verbose:
                msg = cmd
            else:
                msg = 'Running...'
            print("[Find asscoiated ORF]:", msg)
            run_cmd(cmd)
            print("[Find asscoiated ORF]: done!")
            print("="*80)
        
        # step 6: call consensus from reads at sgmRNA split sites
        pos_list = []
        skip_next_step = 0
        f = open(r'{prefix}.{type_prefix}_adjusted_subgenome_pos.tab'.format(prefix=prefix, type_prefix=type_prefix))
        for line in f:
            if not line.startswith('#') and not line.startswith('Adjust'):
                pos_list.append(line.strip().split('\t')[1])
        f.close()
        pos_list = list(set(pos_list))
        if pos_list:
            cmd = '''{py_inter} {script_dir}/call_spliced_site_consensus.py {prefix}.filtered.sam {prefix}.sam.header {pos_list} {prefix} {type_prefix} {py_inter}
            '''.format(py_inter=py_inter, script_dir=script_dir, pos_list=','.join(pos_list), prefix=prefix, type_prefix=type_prefix)
            if verbose:
                msg = cmd
            else:
                msg = 'Running...'
            print("[Call subgenomic mRNA consensus]: ", msg)
            run_cmd(cmd)
            print("[Call subgenomic mRNA consensus]: done!")
            print("="*80)
            # remind the user to provide gff file
            if len(pos_list) <= 4 and not gff:
                print("**Warning**: Only", len(pos_list), "TRS-B positions are found by RNA-seq data. Try to provide an annotation (such as GFF3 format) file to improve the results.")
        else:
            skip_next_step = 1
            print("**Warning**: No TRS-B position is found by RNA-seq data! Only output positions found by potential TRS-B.")
        
        # step 7: format alignment
        # get TRS-L sequence
        if not skip_next_step:
            trs_l_seq = refseq[:150]
            trs_l_file = tmpdir+'/TRS_L.fasta'
            with open(trs_l_file, 'w') as g:
                g.write('>TRS_L\n'+trs_l_seq+'\n')
            sgmrna_doc = tmpdir+'/consensus'
            
            cmd = '''{py_inter} {script_dir}/format_alignment.py {trs_l_file} {sgmrna_doc} {ref_file} > {prefix}.{type_prefix}.format_align.txt
            '''.format(py_inter=py_inter, script_dir=script_dir, trs_l_file=trs_l_file, sgmrna_doc=sgmrna_doc, ref_file=ref, prefix=prefix, type_prefix=type_prefix)
            if verbose:
                msg = cmd
            else:
                msg = 'Running...'
            print("[Format alignment]: ", msg)
            run_cmd(cmd)
            print("[Format alignment]: done!")
            print("="*80)
            
        # step.8: combine only genome and RNA-seq results into a single file
        if gff:
            cmd = '''{py_inter} {script_dir}/combine_only_genome_result.py {prefix}.{type_prefix}_adjusted_subgenome_pos_with_ORF.tab {prefix}.genome_subgenome_pos.tab > {prefix}.{type_prefix}_adjusted_subgenome_pos_with_ORF_combined.tab
            '''.format(py_inter=py_inter, script_dir=script_dir, prefix=prefix, type_prefix=type_prefix)
            if verbose:
                msg = cmd
            else:
                msg = 'Running...'
            print("[Add positions found by sequence similarity]: ", msg)
            run_cmd(cmd)
            print("[Add positions found by sequence similarity]: done!")
            print("="*80)

        # step.9: merge sgmRNA sites and format alignment into output file
        if gff:
            cmd = '''cat {prefix}.TRS_L.txt {prefix}.{type_prefix}_adjusted_subgenome_pos_with_ORF_combined.tab {prefix}.{type_prefix}.format_align.txt > {output}'''.format(output=outputfile, prefix=prefix, type_prefix=type_prefix)
        else:
            cmd = '''cat {prefix}.TRS_L.txt {prefix}.{type_prefix}_adjusted_subgenome_pos.tab {prefix}.{type_prefix}.format_align.txt > {output}'''.format(output=outputfile, prefix=prefix, type_prefix=type_prefix)
        if verbose:
            msg = cmd
        else:
            msg = 'Running...'
        print("[Write output file]: ", msg)
        run_cmd(cmd)
        print("[Write output file]: done!")
        print("="*80)
    
    # step.10: remove tmp files
    cmd = '''rm {prefix}.*'''.format(prefix=prefix)
    print("[Delete temporary files]: ", cmd)
    run_cmd(cmd)
    print("[Delete temporary files]: done!")
    print("="*80)

if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Subgenomic mRNA Position Exploration with RNA-seq (SuPER): Identify TRS-B start positions in reference genome from SAM alignments\n'+
    '**Notice** In order to maximize the speed of SuPER, pypy is recommended as the python interpreter.\n'+
    'Software dependeny: infernal(>=v.1.1.2).\n'+
    'Python modules repuired: biopython, distance, progress.')
    # Essential Input
    parser.add_argument('-r', '--reference', help='reference genome in fasta format', required=True,dest='reference',metavar='')
    parser.add_argument('-g', '--group', help='group of reference genome (e.g. Alpha, Beta, Gamma, Delta)',required=True,dest='group',choices=['Alpha', 'Beta', 'Gamma', 'Delta'])
    parser.add_argument('-a', '--annotate', help='gff3 annnotation for reference genome',required=True,dest='gff',metavar='')

    # Optional Input
    parser.add_argument('-i', '--input', help='input sam file', required=False, dest='input',metavar='')
    parser.add_argument('-d', '--datetype', help='whether sam file data is from Illumina RNA-Seq or Oxford Nanopore data. Can only be 1(for RNA-seq) or 0(for Nanopore Direct RNA Sequencing) [Default:1]',type=int,default=1,required=False,dest='datatype',metavar='')
    parser.add_argument('-l', '--csl', help='the Core Sequence(CS) of leader sequence i.e. TRS-L in reference genome. If not given it will be detected automatically',required=False,dest='csl',metavar='')
    parser.add_argument('-p', '--program', help='the program used to get sam file. If not given it will be detected automatically. **Notice** bwa, minimap2, hisat2 are recommended!', required=False, dest='program')
    parser.add_argument('-c', '--cutoff', help='the cutoff/threshold for the spliced read number so as to keep a spliced site [Default: 3]', type=int, default=3, required=False, dest='cutoff',metavar='')
    parser.add_argument('-v', '--verbose', help='get more information while running in verbose mode', action='store_true', required=False, dest='verbose')
    parser.add_argument('-@', '--threads', help='number of threads used [Default: 1]',type=int, default=1, required=False,dest='core',metavar='')
    # Output
    parser.add_argument('-t', '--tmpdir', help='temporary working directory [Default: tmp]', required=False,dest='tmpdir',default='tmp',metavar='')
    parser.add_argument('-o', '--output', help='output file', required=True, dest='output',metavar='')
    
    args = parser.parse_args()
    
    #check_tools
    for i in ["cmalign","samtools"]:
        if not is_tool(i):
            print("tool {i} is not installed".format(i=i) )
            sys.exit(0)
    
    pipeline(args)
