from sys import argv
from sys import exit
import os
import pandas as pd

from annotation import generate_annotation, process_annotation
from calculation import calculate_score

def getopts(argv):
    opts = {}  
    while argv:  
        if argv[0][0] == '-' and len(argv[0]) > 1 and argv[0][1] != 'h':  
            if len(argv) > 1:
                opts[argv[0]] = argv[1]  
        argv = argv[1:]  
    return opts

if __name__ == '__main__':
    myargs = getopts(argv)
    
    if '-h' in myargs or len(myargs) == 0:  
        print('\nUsage: ./metacompare.py -c filename1.fa [-t 64 -b 1 -s path_prefix -db_dir path] \n')
        print('\t-c: Specify FASTA file containing assembled contigs [required]')
        print('\t-t: Specify the number of threads (default: 64).')
        print('\t-b: Specify pipeline [0: both (default), 1: eco risk, 2: human risk].')
        print('\t-o: Output file path.')
        print('\t-s: Skip annotation. Provide the prefix/path to existing Prodigal files.')
        print('\t-db_dir: Path to the directory metacmpDB (e.g. /path/to/metacmpDB)')
        print()
        exit()

    # --- Setup Database Paths ---
    script_dir = os.path.dirname(os.path.abspath(__file__))
    default_path = os.path.join(script_dir, "metacmpDB")
    
    # metacmp_root is the direct path to the metacmpDB folder
    metacmp_root = myargs.get('-db_dir', default_path)
    
    mge_len_file = os.path.join(metacmp_root, "MGE_len.txt")
    pathogen_file = os.path.join(metacmp_root, "pathogen_list.txt")
    eskape_file = os.path.join(metacmp_root, "eskape.txt")
    gtdb_path = os.path.join(metacmp_root, "GTDB/gtdb")

    nthread = myargs.get('-t', '64')
    out_dir = myargs.get('-o', '')
    pipeline = int(myargs.get('-b', '0'))
    
    sample_name = os.path.splitext(os.path.basename(myargs['-c']))[0]
    out_file = os.path.join(out_dir, sample_name + "_out.txt")    
    
    # --- Annotation Step (Run or Skip) ---
    if '-s' in myargs:
        print(f"Skipping annotation. Using existing files with prefix: {myargs['-s']}")
        prefix = myargs['-s']
        # We fill the list so downstream functions have the file paths they expect
        annotated_data = [
            f"{prefix}.gene.faa", 
            f"{prefix}.gene.fna", 
            f"{prefix}.gff",      
            f"{prefix}.renamed"   
        ]
    else:
        print(f"Running annotation step using database root: {metacmp_root}")
        # Pass the database root so Diamond and MMseqs find their files
        annotated_data = generate_annotation(myargs['-c'], out_dir, nthread, metacmp_root)
    
    # --- Pipeline Execution ---
    if pipeline == 1:        
        data_to_be_processed = [annotated_data[0], annotated_data[1], annotated_data[2]]
        filtered_data = process_annotation(data_to_be_processed, mge_len_file, pathogen_file)
        result = calculate_score(myargs['-c'], filtered_data, pipeline)
    elif pipeline == 2:        
        data_to_be_processed = [annotated_data[3], annotated_data[1], annotated_data[2]]
        filtered_data = process_annotation(data_to_be_processed, mge_len_file, eskape_file)
        result = calculate_score(myargs['-c'], filtered_data, pipeline)
    else:
        # Eco Risk
        data_to_be_processed_e = [annotated_data[0], annotated_data[1], annotated_data[2]]
        filtered_data_e = process_annotation(data_to_be_processed_e, mge_len_file, pathogen_file)
        result_e = calculate_score(myargs['-c'], filtered_data_e, 1)
        
        # Human Risk
        data_to_be_processed_h = [annotated_data[3], annotated_data[1], annotated_data[2]]
        filtered_data_h = process_annotation(data_to_be_processed_h, mge_len_file, eskape_file)
        result_h = calculate_score(myargs['-c'], filtered_data_h, 2)
        
        result = pd.concat([result_e, result_h])
    
    result.to_csv(out_file, header=True, index=None, sep="\t")
    print(f"Done. Results saved to: {out_file}")