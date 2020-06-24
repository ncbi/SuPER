# SuPER
**Su**bgenomic mRNA **P**osition **E**xploration with **R**NA-seq

## Brief Introduction
We developed the tool SuPER to identify TRS-B sites in coronavirus genomes. SuPER first uses a covariance model derived from Rfam to identify TRS-L via profile-based sequence and structure scoring. Then, SuPER identifies TRS-B sites either by identifying template switching junctions with RNA-seq or in the absence of RNA-seq by identifying sequences preceding genes that are similar to the TRS-L CS as putative TRS-B CS.

## Quick Start
```
usage: SuPER.py [-h] -r REFERENCE -g {Alpha,Beta,Gamma,Delta} -a GFF [-i]
                [-d {0,1}] [-l] [-p PROGRAM] [-c] [-v] [-@] [-t] -o

Subgenomic mRNA Position Exploration with RNA-seq (SuPER): Identify TRS-B
start positions in reference genome from SAM alignments **Notice** In order to
maximize the speed of SuPER, pypy is recommended as the python interpreter.
Software dependeny: infernal(>=v.1.1.2). Python modules repuired: biopython,
distance, progress.

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        reference genome in fasta format
  -g {Alpha,Beta,Gamma,Delta}, --group {Alpha,Beta,Gamma,Delta}
                        group of reference genome (e.g. Alpha, Beta, Gamma,
                        Delta)
  -a GFF, --annotate GFF
                        gff3 annnotation for reference genome
  -i , --input          input sam file
  -d {0,1}, --datetype {0,1}
                        whether sam file data is from Illumina RNA-Seq or
                        Oxford Nanopore data. Can only be 1(for RNA-seq) or
                        0(for Nanopore Direct RNA Sequencing) [Default:1]
  -l , --csl            the Core Sequence(CS) of leader sequence i.e. TRS-L in
                        reference genome. If not given it will be detected
                        automatically
  -p PROGRAM, --program PROGRAM
                        the program used to get sam file. If not given it will
                        be detected automatically. **Notice** bwa, minimap2,
                        hisat2 are recommended!
  -c , --cutoff         the cutoff/threshold for the spliced read number so as
                        to keep a spliced site [Default: 3]
  -v, --verbose         get more information while running in verbose mode
  -@ , --threads        number of threads used [Default: 1]
  -t , --tmpdir         temporary working directory [Default: tmp]
  -o , --output         output file
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
##### Position Table
For example:
Adjusted_TRS-B_position|Original_TRS-B_position| Putative_TRS-B|Hamdist_from_TRS-L|Recommend_label*|	Distance_from_ORF|Associated_ORF:start-end|Supported_by_RNAseq|Supported_read_list
---|---|---|---|---|---|---|---|----
20559|20555|AACTAAAT|1|-|11|surface glycoprotein:20570-24091|Yes|ERR3460958.266626.0,ERR3460958.627131.0,...
24468|24468|AACCACAC|2|Not recommended|14|4b protein:24482-24748|No|NA
25672|25665|AACTGAAC|1|-|14|nucleocapsid protein:25686-26855|Yes|ERR3460958.1663.0,ERR3460958.4848.0,...

*Notice: The TRS-B with hamming distance from TRS-L greater than 1 would be labeled as 'Not recommended'.

- Column 1: (most important) TRS-B position adjusted by nearby TRS-B core sequence (Column 3), because sg mRNA might be mapped to the genome with slight shift due to mismatch in alignment algorithm. And TRS-B core sequences are used to adjust the positions. If body core sequence is not found within the original subgenome position surrounding interval of [core seq start position-30bp, core seq end position+30bp], this column is the same as "Original_TRS-B_position".
- Column 2: TRS-B position obtained directly from RNA-Seq data if provided. Otherwise, this position is determined by searching through the reference genome to seek the  TRS-B sequence most similar to TRS-L.
- Column 3: The putative TRS_B core sequence with minimal hamming distance from TRS-L within the original TRS-B position nearby region. 'Not found' means no TRS-B sequence is found.
- Column 4: The hamming distance between TRS-B and TRS-L. 'NA' means no TRS-B is available.
- Column 5: When hamming distance is greater than 2, this record will be labelled as "Not recommended". Yet it is worth noting that not all TRS-B whose hamming distance > 2 should be ignored, and it is recommended to remove them with caution.
- Column 6: The distance from TRS-B start position to the nearest downstream ORF.
- Column 7: If gff annotation file is given, the associated downstream (<170 bp) ORF will be found. Make sure the 'CDS' regions are annotated in the file. 
- Column 8: Whether this position is supported by RNA-seq data. 'No' means this position is found only with a sequence similar to TRS-L.
- Column 9: If this position is supported by RNA-seq data, the spliced reads (shown as 'read_name.sam_flag') are listed.

##### Alignment
After finding TRS-B, the supported reads nearby the position would be collected to reconstruct a consensus subgenomic mRNA, which would be further base-paired with TRS-L an TRS-B in the reference genome.
For example:

```
Position 20555 alignment:
   TRS-L 41 TTTAGACTTTGTGTCTACTTTTCTCAACTAAACGAAATTTTTGCTATG
            ||||||||||||||||||||||||||||||||..||||.||||.|.||
     sgmRNA TTTAGACTTTGTGTCTACTTTTCTCAACTAAATAAAATGTTTGTTTTGCTTGTTGCATATGCCTTGTTGCAT
                                 |||||||||||||||||||||||||||||||||||||||||||||||||||
                     TRS-B 20555 TCTCAACTAAATAAAATGTTTGTTTTGCTTGTTGCATATGCCTTGTTGCAT
```

## Examples

- SuPER can guess the subgenomic mRNA positions even without offering RNA-seq data but with annotation file(i.e. gff3 file):
```
pypy3 SuPER.py -r KY770850.fa -g Alpha -a KY770850.gff3 -t output -o output/KY770850_genome.out.tab
```

- Input Illumina RNA-seq data aligned by bwa and set the supported read number cutoff as 20 (-c 20):
```
pypy3 SuPER.py -i SARS2_bwa_aligned.sam -r MN908947.3.fasta -c 20 -g Beta -a MN908947.3.gff3 -t output -o output/SARS2_bwa_rnaseq.out.tab
```

- Input Illumina RNA-seq data aligned by hisat2:
```
pypy3 SuPER.py -i SARS2_hisat2_aligned.sam -r MN908947.3.fasta -g Beta -a MN908947.3.gff3 -t output -o output/SARS2_hisat2_rnaseq.out.tab
```

- Input Oxford Nanopore RNA-seq data. Tell SuPER the alignment program is 'minimap2' and the data type is Nanopore (-d 0):
```
pypy3 SuPER.py -i SARS2_minimap2_aligned.sam -r MN908947.3.fasta -c 20 -g Beta -p minimap2 -d 0 -a MN908947.3.gff3 -t output -o output/SARS2_minimap2_nanopore.out.tab
```
## Citation
[Yiyan Yang, Wei Yan, Brantley Hall, Xiaofang Jiang. Characterizing transcriptional regulatory sequences in coronaviruses and their role in recombination. bioRxiv. 2020](https://www.biorxiv.org/content/10.1101/2020.06.21.163410v1)