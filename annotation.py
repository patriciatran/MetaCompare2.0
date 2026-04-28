import os, sys, glob
import subprocess
import pandas as pd
import shutil
import warnings

warnings.simplefilter(action='ignore', category=Warning)

IDENTITY = 60
MIN_ALIGN_LENGTH = 25
COVG_OF_ALIGN_MGE = 0.9

def generate_annotation(contig_file, out_dir, nthread='64', metacmp_db_path=None):
    if not out_dir:
        out_dir = './'

    sample_name = os.path.basename(contig_file).split('.')[0]
    ext = os.path.basename(contig_file).split('.')[-1]
    
    prodigal_file = os.path.join(out_dir, f"{sample_name}_genes.{ext}") 
    prodigal_file_faa = os.path.join(out_dir, f"{sample_name}_proteins.{ext}") 

    # 1. Prodigal
    if not os.path.exists(prodigal_file):
        print("Running Prodigal...")
        subprocess.call(["pprodigal", "-i", contig_file, "-d", prodigal_file, "-a", prodigal_file_faa, "-p", "meta", "-C", "4000", "-T", "32", "-o", os.path.join(out_dir, "prodigal_log")])
    else:
        print("Skipping: Prodigal output already exists")

    # 2. Diamond ARGDB
    arg_name = sample_name + "_ARG.csv"        
    if not os.path.exists(os.path.join(out_dir, arg_name)):
        print('Running Diamond Blastx on ARGDB')
        db_path = os.path.join(metacmp_db_path, "ARGDB")
        subprocess.call(["diamond", "blastx", "-d", db_path, "--query", prodigal_file, \
                        "--out", os.path.join(out_dir, arg_name), "--outfmt", "6", \
                        "--threads", nthread, "--evalue", "1e-10"])
    else:
        print('Skipping: Diamond output against ARGs already exists')
        
    # 3. Diamond ARGDB_hh
    arg_name_hh = sample_name + "_hh_ARG.csv"    
    if not os.path.exists(os.path.join(out_dir, arg_name_hh)):
        print('Running Diamond Blastx on ARGDB_hh')
        db_path = os.path.join(metacmp_db_path, "ARGDB_hh")
        subprocess.call(["diamond", "blastx", "-d", db_path, "--query", prodigal_file, \
                        "--out", os.path.join(out_dir, arg_name_hh), "--outfmt", "6", \
                        "--threads", nthread, "--evalue", "1e-10"])
    
    # 4. Diamond MGEDB
    mge_name = sample_name + "_MGE.csv"
    if not os.path.exists(os.path.join(out_dir, mge_name)):
        print('Running Diamond Blastx on MGEDB')
        db_path = os.path.join(metacmp_db_path, "MGEDB")
        subprocess.call(["diamond", "blastx", "-d", db_path, "--query", prodigal_file, \
                        "--out", os.path.join(out_dir, mge_name), "--outfmt", "6", \
                        "--threads", nthread, "--evalue", "1e-10"])

    # 5. MMSeqs2 (GTDB)
    pathogen_name = sample_name + "_Pathogens.tsv"
    if not os.path.exists(os.path.join(out_dir, pathogen_name)):
        print('Running mmseq2 on GTDB')
        gtdb_full_path = os.path.join(metacmp_db_path, "GTDB/gtdb")
        
        # Format out_dir for the shell script
        sh_out_dir = out_dir if out_dir.endswith('/') else out_dir + "/"
        sh_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "mmseq.sh")
        
        subprocess.call(['sh', sh_script, contig_file, sh_out_dir, pathogen_name, gtdb_full_path])
        
        try:
            shutil.rmtree(os.path.join(out_dir, "sample.tmpFolder"), ignore_errors=True)
            for f in glob.glob(os.path.join(out_dir, "sample.contigs*")):
                os.remove(f)
            for f in glob.glob(os.path.join(out_dir, "sample.assignments*")):
                os.remove(f)
        except OSError as e:
            print("Error cleaning up MMSeqs files: %s" % e.strerror)
    
    return [os.path.join(out_dir, arg_name), os.path.join(out_dir, mge_name), os.path.join(out_dir, pathogen_name), os.path.join(out_dir, arg_name_hh)]

def filter_diamond(filename):
    data = pd.read_csv(filename, sep='\t', header=None)
    data.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
    filtered_data = data[(data.identity > IDENTITY) & (data.alignLen > MIN_ALIGN_LENGTH)]
    return filtered_data

def process_annotation(data_files, mge_len_file, pathogen_list):    
    arg_file, mge_file, path_file = data_files[0], data_files[1], data_files[2]
    
    # Process ARGs
    arg_data = filter_diamond(arg_file) if os.path.getsize(arg_file) > 0 else pd.DataFrame()
        
    # Process MGEs
    if os.path.getsize(mge_file) > 0:
        mge_data = filter_diamond(mge_file) 
        mge_len = pd.read_csv(mge_len_file, sep='\t', header=None, names=['sub_id', 'ref_gene_leng'])
        mge_merged = pd.merge(mge_data, mge_len, on='sub_id', how='left')
        mge_final = mge_merged[mge_merged.alignLen > (mge_merged.ref_gene_leng * COVG_OF_ALIGN_MGE)]        
    else:
        mge_final = pd.DataFrame()
            
    # Process Pathogens
    if os.path.getsize(path_file) > 0:
        path_data = pd.read_csv(path_file, sep='\t', header=None)
        path_data.columns = ['id', 'NCBI_ID','rank', 'name', 'nPass', 'nRetain', 'nAssign', 'bit', 'taxonomy']
        path_filtered = path_data[path_data['rank'].isin(["family", "genus", "species", "strain"])]
        
        pathogens_ref = pd.read_csv(pathogen_list, sep="\t")
        ranks = pathogens_ref['rank'].unique()
        path_final_list = []

        for rank in ranks:
            valid_names = pathogens_ref[pathogens_ref['rank'] == rank]["name"].tolist()
            match = path_filtered[(path_filtered['rank'] == rank) & (path_filtered['name'].isin(valid_names))]
            path_final_list.append(match)
        
        # Special strain handling
        strains = path_filtered[path_filtered['rank'] == "strain"].copy()
        if not strains.empty:
            strains['short_name'] = strains['name'].apply(lambda x: ' '.join(x.split()[:2]))
            valid_species = pathogens_ref[pathogens_ref['rank'] == "species"]["name"].tolist()
            match_strains = strains[strains['short_name'].isin(valid_species)].drop(columns=['short_name'])
            path_final_list.append(match_strains)
            
        path_all = pd.concat(path_final_list).drop_duplicates()
    else:
        path_all = pd.DataFrame()
        
    return [arg_data, mge_final, path_all]