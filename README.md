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
