Directory structure:
for downloading data and parsing gtf
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


Directory structure for Motif discovery in eCLIP using PRIESSTESS + STREME + Tomtom
```
eclip-priesstess-motif-pipeline/
├─ README.md
├─ LICENSE
├─ .gitignore
├─ environment/
│  ├─ conda.yml
│  └─ requirements.txt
├─ config/
│  ├─ config.yaml
│  └─ samples.tsv
├─ data/
│  ├─ raw/
│  │  └─ eclip/                        # raw eCLIP peaks/merged tables (gz/csv/bed)
│  ├─ processed/
│  │  ├─ pos_bed/                      # positive peak beds
│  │  ├─ neg_bed/                      # negative/control beds
│  │  ├─ fasta/                        # sequences extracted from bed (pos/neg)
│  │ 
│  └─ references/
│     ├─ genome/                       # hg38 fasta + index files (not committed)
│     └─ motifs/
│        ├─ RBPDB.meme                 # motif library for Tomtom
│        └─ JASPAR.meme                # optional (if used)
├─ scripts/
│  ├─ 01_preprocess_split_pos_neg.sh
│  ├─ 02_bed_to_fasta.sh
│  ├─ 03_run_priesstess.sh
│  ├─ 04_extract_streme_motifs.py
│  ├─ 05_streme_to_meme.py
│  ├─ 06_run_tomtom.sh
│  └─ 07_summarize_tomtom.py
├─ notebooks/
│  └─ explore_results.ipynb
├─ results/
│  ├─ priesstess/
│  ├─ streme/
│  ├─ tomtom/
│  └─ summary/
│     ├─ combined_tomtom.tsv
│
└─ logs/
```

  

## Overview
This repository provides a step-by-step guide for processing eCLIP (enhanced CrossLinking and ImmunoPrecipitation) data from the ENCODE database and running the PRIESSTESS + STREME + Tomtom models for RNA-binding protein (RBP) motif identification.

## Requirements

RNAfold (PRIESSTESS was developed with version 2.4.11)

STREME (PRIESSTESS was developed with version 5.3.0)

python3 (PRIESSTESS was developed with version 3.8)

sklearn (PRIESSTESS was developed with version 0.23.2)

skopt python package (PRIESSTESS was developed with version 0.8.1)

**About the dataset** :
This is a eCLIP (enhanced CrossLinking and ImmunoPrecipitation) data obtained from the ENCODE (Encyclopedia of DNA Elements) database. eCLIP is a high-throughput technique used to map RNA-binding protein (RBP) interaction sites across the transcriptome. The current dataset focuses on K562 and HepG2 celllines. 

**File formats** :
FASTQ (.fastq.gz): Raw sequencing reads before processing   
BED (.bed): Genomic coordinates of peaks (binding sites)      

## eCLIP data preprocessing
This section describes the preprocessing steps used to generate sequence-level inputs from BED files for downstream analysis. All commands were executed in a Unix/Linux environment.
   - Start from a raw eCLIP peak file/table (BED/CSV).
   - Generate **positive** (enriched peaks) and **negative** (background/control) sets.
   - Convert the bed files into positive and negative sets where the positive set contains high-scoring regions. The negative set has similar regions but shifted away, which can be used as a control.

    1. Extract Sequences from BED Files
    
    Genomic sequences were extracted from the reference genome using BEDTools based on genomic coordinates provided in BED files.
    
    ```bash
    bedtools getfasta -fi genome.fa -bed input.bed -fo output.fa
    ```
    -fi: Reference genome FASTA file
    -bed: Input BED file containing genomic intervals
    -fo: Output FASTA file containing extracted sequences
    
    2. Remove FASTA Headers
    
    To retain only the raw nucleotide sequences (without chromosome names or coordinates), FASTA headers were removed using awk.
    ```bash
    awk '!/^>/' ENCFF002HYO_negative.fa > ENCFF002HYO_neg.fa
    ```
    This step produces a plain sequence file with one sequence per line.
    
    3. Convert BED to TXT Without Headers (Direct Extraction)
    
    In some cases, sequences were extracted directly from BED files and written as header-free text files.
    ```bash
    bedtools getfasta -fi GRCh38.primary_assembly.genome.fa -bed ENCFF227EJF_positive.bed | grep -v "^>" > ENCFF227EJF_positive.txt
    ```
    This approach combines sequence extraction and header removal in a single step.
    
    4. Delete Locked Directories
    
    During preprocessing, certain output directories were locked due to active .nfs* files. These were identified and resolved before deletion.
    ```bash
    lsof PRIESSTESS_output/.nfs2754c561ee16cb910000cee2
    kill -9 <PID>
    ```
    Once the locking process was terminated, the directory could be safely removed.
    
    5. Convert DNA Sequences to RNA
    
    DNA sequences were converted to RNA by replacing thymine (T) with uracil (U).
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
**Generated Files**

ENCFF002HYO_neg.fa: Negative sequences without headers.

ENCFF227EJF_positive.txt: Positive sequences extracted from BED files.

ENCFF031FMO_rna_positive.txt: RNA sequences after T→U conversion.

## Run PRIESSTESS
   - Execute PRIESSTESS on pos/neg inputs.
   - Produce **STREME** outputs (motifs) and PRIESSTESS scan results.
 
```
PRIESSTESS -fg foreground_file -bg background_file
```
/work/talisman/smuthyala/motif_identification/PRIESSTESS_for_eClip/PRIESSTESS/PRIESSTESS -fg /work/talisman/smuthyala/motif_identification/PRIESSTESS_for_eClip/eCLIP-ENCODE-GRCH38/data/eCLIP/ENCFF031FMO_rna_positive.txt -bg /work/talisman/smuthyala/motif_identification/PRIESSTESS_for_eClip/eCLIP-ENCODE-GRCH38/data/eCLIP/ENCFF031FMO_rna_negative.txt

## PRIESSTESS OUTPUT 

1. Model & Training Files
PRIESSTESS_model.sav → The trained PRIESSTESS model (serialized using Python's pickle or similar).
PRIESSTESS_model_weights.tab → The learned weights for sequence-structure features used in classification.
PRIESSTESS_arguments.txt → Stores the command-line arguments used to run PRIESSTESS (input/output files, parameters, etc.).
LR_training_set.tab → The training dataset used for logistic regression (likely contains feature vectors for RNA sequences).

2. Feature Representation & Annotation Files
annotation_alphabets_header.tab → Defines the annotation format for sequence and structure representations.
bg_alphabet_annotations.tab → Alphabet annotations for background (negative) RNA sequences.
fg_alphabet_annotations.tab → Alphabet annotations for foreground (positive) RNA sequences.

3. Sequence-Structure Feature Files
These files contain numerical encodings of RNA sequences combined with structural information at different resolutions:
seq-4, seq-struct-8, seq-struct-16, seq-struct-28 → Encoded sequence-structure features at different resolutions.
struct-2, struct-4, struct-7 → Encoded structure-only features at different resolutions.

* What are these seq-struct files? PRIESSTESS extracts RNA sequence-structure features using different context lengths (4, 8, 16, 28, etc.). Larger numbers indicate a broader window of sequence and structural context considered for learning.

4. Evaluation & Performance Metrics
heldout_data.tab → Contains test/validation data that was held out from training.
test_PRIESSTESS_model_ON_heldout_auroc.tab → Contains the AUROC (Area Under the Receiver Operating Characteristic curve) for evaluating the trained model's performance.##

### 1. Sequence Data Files
fg_LR.fa → Foreground (positive) sequences used for Logistic Regression (likely FASTA format).

bg_LR.fa → Background (negative) sequences used for Logistic Regression.

fg_test.fa → Foreground test sequences (held-out data for validation/testing).

bg_test.fa → Background test sequences (held-out data for validation/testing).

🔍 Interpretation:
These files contain RNA sequences in FASTA format, split into foreground (binding sites) and background (non-binding sites) for training/testing PRIESSTESS.

### 2. Motif Discovery Files
fg_STREME.fa → Foreground sequences formatted for STREME (motif discovery tool).

bg_STREME.fa → Background sequences formatted for STREME.

streme.html → HTML report with visualization of discovered motifs.

streme.txt → Text file summarizing motif discovery results.

streme.xml → XML-formatted motif results.

### Interpretation:
These files come from STREME, a motif discovery algorithm.

Foreground motifs (fg_STREME.fa): Represent common sequence patterns in bound RNA sites.

Background motifs (bg_STREME.fa): Used to distinguish real motifs from random patterns.

Results (streme.*): Contain significant sequence motifs found in the foreground dataset.

### Key Takeaways
The FASTA files store RNA sequences used for classification and motif discovery.

The STREME output files provide insights into sequence motifs enriched in RNA-protein binding sites.

The motif discovery results can help interpret RNA-binding protein preferences.

## **Extract motifs from PRIESSTESS output**

Goal: standardize motifs into a format downstream tools (like Tomtom) can read and curate a high-confidence set.

- Parse streme.xml/streme.txt.
- Convert to **MEME** motif format

**Outputs**

- A curated discovered.meme containing all retained PWMs (+ optional structure metadata).
- A companion table with motif ID, width, E-value, #sites, info content.

## **Run Tomtom against known RBP motif databases**

Goal: map your discovered motifs to known RBP binding models to propose candidate RBPs.

Following motif discovery with **RBPDB (RNA-Binding Protein DataBase)**, we automated the mapping of candidate motifs to known RNA-binding proteins:

1. We parse the output motif files and query RBPDB's.
2. The script batch-submitted motifs to RBPDB's backend and retrieved RBPs with known affinity for similar sequences based on similarity scoring.
3. We compiled the matches into an annotation table linking each discovered motif to potential RBP identities, including confidence scores and alignment metrics.

**Outputs**

- tomtom.tsv with: query_motif, target_motif (RBP), p-value, q-value, overlap, offset, orientation.
- Ranked matches per motif; often several plausible RBPs per motif.
  
## Notes

1. Ensure that the reference genome (genome.fa or GRCh38.primary_assembly.genome.fa) is correctly indexed before running bedtools getfasta.
2. The kill -9 command should be used cautiously to terminate processes.
3. Confirm sequence integrity after sed modifications before proceeding with PRIESSTESS.

