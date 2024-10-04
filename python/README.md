```{sh}
pip -m venv .eclip
source .eclip/bin/activate
pip install -r requirements.txt
python eCLIP_parser.py
```

If you want to chage the GTF version please do so in line 75 of ```eCLIP_parser.py```. It is set to v29 in the the code.
