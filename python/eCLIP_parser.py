#!/usr/bin/env python
# coding: utf-8

from gtfparse import read_gtf
import polars as pl
import glob
import urllib.request as request
import os
import logging
import sys
import git

def trim_chr(x):
    ''' function for aligning the chromosomes in the GTF with the eCLIP data
        arguments:
            x: chromosome name
        return:
            formatted chromosome name
    '''
    sp = x.split('_')
    return sp[0] if len(sp) == 1 else sp[1].replace('v', '.')


def get_git_root(path):
    ''' gets the root of the git repository
    '''
    git_repo = git.Repo(path, search_parent_directories=True)
    git_root = git_repo.git.rev_parse("--show-toplevel")
    return git_root
        

def init_logging(file):
    ''' initialize the logger
        argument:
            file: name of the file to log the output to
    '''
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
        
    logging.basicConfig(filename=file, filemode='a', format='%(levelname)s [%(asctime)s]: %(message)s', level=logging.INFO, datefmt='%I:%M:%S')
    st = logging.StreamHandler(sys.stdout)
    st.setFormatter(logging.Formatter('%(levelname)s [%(asctime)s]: %(message)s', datefmt='%I:%M:%S'))
    logging.getLogger().addHandler(st)
        
        
def main():
    
    verbose = False
    
    init_logging('eCLIP_ENCODE.log')
    
    git_path = get_git_root('.')
    
    # read the metadata file
    logging.info(f'reading metadata file from {git_path}/data/eCLIP/metadata.tsv for list of eCLIP experiments')
    df_m = pl.read_csv(git_path+'/data/eCLIP/metadata.tsv', separator='\t')
    file_list = df_m.filter((pl.col('File assembly') == 'GRCh38') & (pl.col('Biological replicate(s)') == '1, 2')).select(pl.col('File download URL')).to_numpy().flatten()

    # download the eCLIP files
    logging.info(f'downloading eCLIP bed files from ENCODE and saving to {git_path}/data/eCLIP/')
    i = 0
    for file in file_list:
        try:
            path = git_path+'/data/eCLIP/'+file.split('/')[-1]
            if not os.path.isfile(path):
                request.urlretrieve(file, path)
                i += 1
                if verbose:
                    logging.info('Downloaded file no. {}: {}'.format(i, file.split('/')[-1]))
        except Exception as e:
            logging.error('Cannot download file {} due to {}'.format(file.split('/')[-1], e))


    # download the correct gtf version
    gtf_version = '29'
    logging.info(f'downloading primary assembly GTF v{gtf_version} to {git_path}/data/gtf')
    gtf_url = 'https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_'+gtf_version+'/gencode.v'+gtf_version+'.primary_assembly.annotation.gtf.gz'
    gtf_path = git_path+'/data/gtf/gencode.v'+gtf_version+'.primary_assembly.annotation.gtf'
    if not os.path.isfile(gtf_path):
        request.urlretrieve(gtf_url, gtf_path)
    df_f = read_gtf(gtf_path)
    df_f = df_f.filter(pl.col('feature') == 'transcript').select(pl.col('seqname', 'start', 'end', 'strand', 'gene_id', 'transcript_id'))


    # read all the eCLIP bed files
    logging.info(f'reading all the bed files')
    queries = []
    l_cols = ["chr", "start", "stop", "dataset_label", "1000", "strand", 
               "log2(eCLIP fold-enrichment over size-matched input)", 
               "-log10(eCLIP vs size-matched input p-value)", "-1", "-1.1"]
    for file in glob.glob("../data/eCLIP/*.bed.gz"):
        q = pl.read_csv(file, has_header=False, separator='\t')
        q.columns = l_cols
        f_acc = file.split('/')[-1].split('.')[0]
        c_line= df_m.filter(pl.col('File accession') == f_acc).select(pl.col('Biosample term name')).to_numpy()[0][0]
        rbp = df_m.filter(pl.col('File accession') == f_acc).select(pl.col('Experiment target')).to_numpy()[0][0].split('-')[0]
        q = q.with_columns(RBP = pl.lit(rbp))
        q = q.with_columns(cell_line = pl.lit(c_line))
        # fix the "." problem in 27 experiments
        if q.filter(pl.col('dataset_label') == '.').shape[0] > 0:
            q = q.with_columns(pl.lit(rbp+'_'+c_line+'_.').alias('dataset_label'))
        queries.append(q)
    df_e = pl.concat(queries)
    df_e = df_e.with_columns(pl.col('chr').map_elements(trim_chr, return_dtype=pl.String))


    # compare the chromosomes in the eCLIP data with that in the GTF
    compare = False
    if compare:
        l_chr_e = df_e.select(pl.col('chr')).unique().to_numpy().flatten()
        l_chr_f = df_f.select(pl.col('seqname')).unique().to_numpy().flatten()
        if verbose:
            logging.info('chr in eCLIP but not in GTF: {}'.format(set(l_chr_e) - set(l_chr_f)))
            logging.info('chr in GTF but not in eCLIP: {}'.format(set(l_chr_f) - set(l_chr_e)))

    # blend in the eCLIP data with the GTF file
    logging.info(f'merging eCLIP data with GTF annotations')
    queries = []
    for row in df_f.iter_rows(named=True):
        df_tmp = df_e.filter((pl.col('start') >= row['start']) 
                             & (pl.col('stop') <= row['end']) 
                             & (pl.col('chr') == row['seqname']) 
                             & (pl.col('strand') == row['strand']))
        if df_tmp.shape[0] != 0:
            df_tmp = df_tmp.with_columns(featureStart = pl.lit(row['start']))
            df_tmp = df_tmp.with_columns(featureEnd = pl.lit(row['end']))
            df_tmp = df_tmp.with_columns(frame = pl.lit(row['strand']))
            df_tmp = df_tmp.with_columns(ENSG = pl.lit(row['gene_id'].split('.')[0]))
            df_tmp = df_tmp.with_columns(ENST = pl.lit(row['transcript_id'].split('.')[0]))
            queries.append(df_tmp)

    df_a = pl.concat(queries)

    # save results in a zipped csv
    logging.info(f'saving annotated eCLIp file to: {git_path}/data/eCLIP/eCLIP_ENCODE_merged_April_2024_GRCh38_GENCODEv{gtf_version}.csv.gz')
    df_a.to_pandas().to_csv(git_path+'/data/eCLIP/eCLIP_ENCODE_merged_April_2024_GRCh38_GENCODEv'+gtf_version+'.csv.gz', index=False)

    
if __name__ == "__main__":
    # execute only if run as a script
    main()
