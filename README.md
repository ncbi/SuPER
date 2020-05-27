# SuPER
## Subgenome Position Exploration by RNA sequencing data

## Quick Start
```
usage: SuPER.py [-h] -r  -g  [-i] [-d] [-a] [-l] [-p] [-c] [-@] [-t] -o

Subgenome Position Exploration by RNA sequencing data (SuPER): Identify
subgenome start positions in reference genome from SAM alignments **Notice**
In order to maximize the speed of SuPER, pypy is recommended as the python
interpreter. Software dependeny: infernal(>=v.1.1.2). Python modules repuired:
biopython, distance, progress.

optional arguments:
  -h , --help         show this help message and exit
  -r , --reference   reference genome in fasta format
  -g , --group       group of reference genome, e.g. Alpha, Beta, Gamma, Delta
  -i , --input       input sam file
  -d , --datetype    whether sam file data is from Illumina RNA-Seq or Oxford
                     Nanopore data. Can only be 1(for RNA-seq) or 0(for
                     Nanopore Direct RNA Sequencing) [Default:1]
  -a , --annotate    gff3 annnotation for reference genome
  -l , --csl         the Core Sequence(CS) of leader sequence i.e. TRS-L in
                     reference genome. If not given it will be detected
                     automatically
  -p , --program     the program used to get sam file. If not given it will be
                     detected automatically. **Notice** bwa, minimap2, hisat2
                     are recommended!
  -c , --cutoff      the cutoff/threshold for the spliced read number so as to
                     keep a spliced site [Default: 3]
  -@ , --threads     number of threads used [Default: 1]
  -t , --tmpdir      temporary working directory [Default: tmp]
  -o , --output      output file
```

## Environment
Although regular Python Interpreter can run SuPER, to use pypy3 as Python Interpreter can speed up SuPER.   
Make sure you install Anaconda first, and you can install pypy3 as follows:

```
conda create -n py3 python=3.6
conda activate py3
conda install -c conda-forge pypy3.6
# make pip available in pypy3
pypy3 -m ensurepip
```
To install dependency like Infernal:

```
conda install -c bioconda infernal=1.1.2
```
To install python modules:

```
pip install biopython==1.70 distance progress
```

### Input
We recommend you use Illumina RNA-Seq data as NGS data and Direct RNA Sequencing (DRS) of Oxford Nanopore as 3GS data to get the sam file.
For metatranscriptomic data, we recommend that host contamination be removed first before getting the final SAM file in order to reduce run time and increase accuracy.
The group classification (i.e. Alpha, Beta, Gamma or Delta) of a coronavirus strain is required.

### Output
Adjusted_subgenome_position | Original_subgenome_position | Putative_body_CS | Hamdist_from_leader_CS |Associated_ORF:start-end | Supported_read_list
---|---|---|---|---|----
25380 |	25382 | AACGAAC | 2 | gene-ORF3a:25393-26220 | SRR11059944.sra.749567.153,SRR11059944.sra.750829.163,...

- Column 1: (most important) Subgenome position adjusted by near body core sequence, because sg mRNA might be mapped to the genome with slight shift due to mismatch in alignment algorithm. And body core sequences are used to adjust the positions. If body core sequence is not found within the original subgenome position surrounding interval of [core seq start position-10bp, core seq end position+15bp], then this column is the same as Original_subgenome_position.
- Column 2: Subgenome position obtained directly from RNA-Seq data.
- Column 3: The body core sequence is found within the original subgenome position surrounding region (<100 bp). 'Not found' means no body core sequence is found.
- Column 4: The hamming distance between the body core sequence and leader core sequence. 'NA' means no body core sequence is found.
- Column 5: If gff annotation file is given, SuPER will find the associated downstream (<100bp) ORF of adjusted_subgenome_position. Make sure 'CDS' regions are annotated in this file. 
- Column 6: The spliced reads (shown as 'read_name.sam_flag') that are found being spliced at Original_subgenome_position.

## Examples

- SuPER can guess the subgenome positions even without offering RNA-seq data but with annotation file(i.e. gff3 file):
```
pypy3 SuPER.py -r KY770850.fa -g Alpha -a KY770850.gff3 -t output -o output/KY770850_genome.out.tab
```

- Even without annotation file:
```
pypy3 SuPER.py -r KY770850.fa -g Alpha -t output -o output/KY770850_genome.out.tab
```

- Input Illumina RNA-seq data aligned by bwa and set the read number cutoff as 20 (-c 20):
```
pypy3 SuPER.py -i SARS2_bwa_aligned.sam -r MN908947.3.fasta -c 20 -g Beta -a MN908947.3.gff3 -t output -o output/SARS2_bwa_rnaseq.out.tab
```

- Input Oxford Nanopore RNA-seq data. Tell SuPER the alignment program is 'minimap2' and the data type is Nanopore (-d 0):
```
pypy3 SuPER.py -i SARS2_bwa_aligned.sam -r MN908947.3.fasta -c 20 -g Beta -p minimap2 -d 0 -a MN908947.3.gff3 -t output -o output/SARS2_minimap2_nanopore.out.tab
```
