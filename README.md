Directory structure:

```
.
├── data
│   ├── eCLIP
│   │   ├── ENCFF002HYO.bed.gz
│   │   ├── ENCFF005GER.bed.gz
│   │   ├── ENCFF007VPS.bed.gz
│   │   ├── ENCFF011NAP.bed.gz
.   .   .
.   .   .
.   .   .
│   │   └── metadata.tsv
│   ├── eCLIP-all
│   │   ├── ENCODE_files.txt
│   │   └── geteCLIPfromENCODE.sh
│   └── gtf
│       └── gencode.v29.primary_assembly.annotation.gtf
├── LICENSE
└── python
    ├── eCLIP_parser.ipynb
    ├── eCLIP_parser.py
    └── requirements.txt
```

__data/eCLIP__: contains all the eCLIP file from ENCODE  
__data/eCLIP-all__ is a placeholder for pull all replicates and merged files for the eCLIP experiments  
__python__: code for downloading and parsing gtfs and eCLIP data  

## Overview
This repository provides a step-by-step guide for processing eCLIP (enhanced CrossLinking and ImmunoPrecipitation) data from the ENCODE database and running the PRIESSTESS model for RNA-binding protein (RBP) motif identification.

## Requirements

RNAfold (PRIESSTESS was developed with version 2.4.11)

STREME (PRIESSTESS was developed with version 5.3.0)

python3 (PRIESSTESS was developed with version 3.8)

sklearn (PRIESSTESS was developed with version 0.23.2)

skopt python package (PRIESSTESS was developed with version 0.8.1)

## About the dataset
This is a eCLIP (enhanced CrossLinking and ImmunoPrecipitation) data obtained from the ENCODE (Encyclopedia of DNA Elements) database. eCLIP is a high-throughput technique used to map RNA-binding protein (RBP) interaction sites across the transcriptome. The current dataset focuses on K562 and HepG2 celllines. 

## File formats 
FASTQ (.fastq.gz): Raw sequencing reads before processing   
BED (.bed): Genomic coordinates of peaks (binding sites)      

## Data Pre processing 
Convert the bed files into positive and negative sets where the positive set contains high-scoring regions. The negative set has similar regions but shifted away, which can be used as a control.
Script for it is it bash scripts folder 

1. Extract Sequences from BED File

Extract sequences from the reference genome using BEDTools:

```bash
bedtools getfasta -fi genome.fa -bed input.bed -fo output.fa
```
2. Remove Chromosome Numbers and Coordinates

To retain only the sequences in the FASTA file without headers:
```bash
awk '!/^>/' ENCFF002HYO_negative.fa > ENCFF002HYO_neg.fa
```
3. Convert BED to TXT Without Headers

Extract sequences and remove headers:
```bash
bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed ENCFF227EJF_positive.bed | grep -v "^>" > ENCFF227EJF_positive.txt
```
4. Delete a Directory

Check and terminate processes locking a directory before deletion:
```bash
lsof PRIESSTESS_output/.nfs2754c561ee16cb910000cee2
kill -9 <PID>
```
5. Convert DNA Sequences to RNA

Replace thymine (T) with uracil (U) in sequences:
```bash
sed 's/T/U/g' dna_sequences.txt > rna_sequences.txt
```
For specific files:
```bash
sed 's/T/U/g' ENCFF031FMO_positive.txt > ENCFF031FMO_rna_positive.txt
```
6. Count the Number of Sequences in Files

Check the number of sequences in multiple files within a directory:
```bash
wc -l *_positive_rna.txt
```
## Processed Files

ENCFF002HYO_neg.fa: Negative sequences without headers.

ENCFF227EJF_positive.txt: Positive sequences extracted from BED files.

ENCFF031FMO_rna_positive.txt: RNA sequences after T→U conversion.

# running PRIESSTESS model 
```
PRIESSTESS -fg foreground_file -bg background_file
```
/work/talisman/smuthyala/motif_identification/PRIESSTESS_for_eClip/PRIESSTESS/PRIESSTESS -fg /work/talisman/smuthyala/motif_identification/PRIESSTESS_for_eClip/eCLIP-ENCODE-GRCH38/data/eCLIP/ENCFF031FMO_rna_positive.txt -bg /work/talisman/smuthyala/motif_identification/PRIESSTESS_for_eClip/eCLIP-ENCODE-GRCH38/data/eCLIP/ENCFF031FMO_rna_negative.txt

# PRIESSTESS Output structure 
PRIESSTESS Output Package

Contents:
- motifs/: Contains discovered motifs in MEME format.
- performance/: AUC, ROC, and accuracy metrics for model evaluation.
- sequences/: Binding site predictions in BED/FASTA format.
- visualizations/: Motif logos and other plots.
- model/: Trained model weights.
- logs/: Run logs and parameters.

## Notes

Ensure that the reference genome (genome.fa or GRCh38.primary_assembly.genome.fa) is correctly indexed before running bedtools getfasta.

The kill -9 command should be used cautiously to terminate processes.

Confirm sequence integrity after sed modifications before proceeding with PRIESSTESS.


ls -d */ | wc -l
