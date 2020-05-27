# SuPER
## Subgenome Position Exploration by RNA sequencing data

> usage: SuPER.py [-h] -i  -r  [-d] [-g] [-l] [-t] -o
> 
> Subgenome Position Exploration by RNA sequencing data (SuPER): Identify
> subgenome start positions in reference genome from SAM alignments 
> **Notice** In order to maximize the speed of SuPER, pypy is recommended as the python
> interpreter. 
> Python modules repuired: biopython, distance.
> 
> optional arguments:  
>   -h, --help         show this help message and exit  
>   -i , --input       input sam file  
>   -r , --reference   reference genome in fasta format  
>   -d , --datetype    if sam file data is from Illumina RNA-Seq or Oxford
>                      Nanopore data. Can only be 1(for RNA) or 0(for Nanopore)  
>   -g , --gff         gff3 annnotation for reference genome  
>   -l , --csl         the Core Sequence(CS) of leader sequence in reference
>                      genome. If not given it will be inferred by program  
>   -t , --tmpdir      tmp dir  
>   -o , --output      output file

**Notice**: 
Use pypy3 as Python Interpreter can speed up SuPER.   
We recommend you use Illumina RNA-Seq data as NGS data and Direct RNA Sequencing (DRS) of Oxford Nanopore as 3GS data.  
For metatranscriptomic data, it will be better to remove host contamination first before getting the final SAM file in order to reduce run time.

### Output
Adjusted_subgenome_position | Original_subgenome_position | Putative_body_CS | Hamdist_from_leader_CS | Supported_read_list | Associated_ORF:start-end
---|---|---|---|---|----
25380 |	25382 | CATAAAC(GAAC) | 2 | SRR11059944.sra.749567.153,SRR11059944.sra.750829.163,... | gene-ORF3a:25393-26220

- Column 1: (most important) Subgenome position adjusted by near body core sequence, because sg mRNA might be mapped to the genome with slight shift due to mismatch in alignment algorithm. And body core sequences are used to adjust the positions. If body core sequence is not found within the original subgenome position surrounding interval of [core seq start position-10bp, core seq end position+15bp], then this column is the same as Original_subgenome_position.
- Column 2: Subgenome position obtained directly from RNA-Seq data.
- Column 3: The body core sequence (often followed by 4-mer flank sequence--'GAAC') is found within the original subgenome position surround region. 'Not found' means no body core sequence is found.
- Column 4: The hamming distance between the body core sequence and leader core sequence. 'NA' means no body core sequence is found.
- Column 5: The spliced reads (shown as 'read_name.sam_flag') that are found being spliced at Original_subgenome_position.
- Column 6: If gff annotation file is given, SuPER will find the associated downstream (<100bp) ORF of adjusted_subgenome_position. Make sure 'CDS' regions are annotated in this file. 

### Examples:

```
# use pypy3 to test Illumina RNA-Seq data with gff file
WD=/home/qiime/RNA-seq_illumina
pypy3 SuPER.py -i $WD/hisat2_output.sam -r $WD/MN908947.3.fna -g $WD/MN908947.3.gff3 -t test_rnaseq -o final.rnaseq.tab
```

```
# use python to test Illumina RNA-Seq data without gff file
WD=/home/qiime/RNA-seq_illumina
python SuPER.py -i $WD/hisat2_output.sam -r $WD/MN908947.3.fna -t test_rnaseq_wo_gff -o final_wo_gff.rnaseq.tab
```

```
# use pypy3 to test Oxford Nanopore data with gff file
WD=/home/qiime/nanopore_seq
pypy3 SuPER.py -i $WD/minimap2_output.sam -r $WD/MT007544.1.fasta -g $WD/MT007544.1.gff3 -d 0 -t test_nanopore -o final.nanopore.tab
```

```
# use python to test Oxford Nanopore data without gff file 
WD=/home/qiime/nanopore_seq
python SuPER.py -i $WD/minimap2_output.sam -r $WD/MT007544.1.fasta -d 0 -t test_nanopore_wo_gff -o final_wo_gff.nanopore.tab
```
