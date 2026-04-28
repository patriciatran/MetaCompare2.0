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
        print('\nUsage: ./metacompare.py -c filename1.fa [-t 64 -b 1 -s path_prefix] \n')
        print('\t-c: Specify FASTA file containing assembled contigs [required]')
        print('\t-t: Specify the number of threads will be used in executing blast (default: 64).')
        print('\t-b: Specify the pipeline to execute [0: both (default), 1: ecological risk score, 2: human health risk score ].')
        print('\t-o: Output file path.')
        print('\t-s: Skip annotation. Provide the prefix/path to existing Prodigal files.')
        print()
        exit()

    # Default values
    nthread = myargs.get('-t', '64')
    out_dir = myargs.get('-o', '')
    pipeline = int(myargs.get('-b', '0'))
    
    mge_len_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/MGE_len.txt")
    sample_name = os.path.splitext(os.path.basename(myargs['-c']))[0]
    out_file = os.path.join(out_dir, sample_name + "_out.txt")    
    
    # Logic to skip or run annotation
    if '-s' in myargs:
        print(f"Skipping annotation. Using existing files with prefix: {myargs['-s']}")
        # We simulate the return list of generate_annotation. 
        # Typically these are paths to .faa, .fna, .gff, etc.
        # Ensure these indices match what your generate_annotation function usually returns.
        prefix = myargs['-s']
        annotated_data = [
            f"{prefix}.gene.faa", # Index 0: proteins
            f"{prefix}.gene.fna", # Index 1: nucleotide genes
            f"{prefix}.gff",      # Index 2: coords
            f"{prefix}.renamed"   # Index 3: specific for pipeline 2 if applicable
        ]
    else:
        print("Running annotation step...")
        annotated_data = generate_annotation(myargs['-c'], out_dir, nthread)    
    
    # Pipeline execution
	# Option 1 is ecological risk only, Option 2 is human risk score only. Otherwise, calculate both scores.
    if pipeline == 1:        
        data_to_be_processed = [annotated_data[0], annotated_data[1], annotated_data[2]]
        pathogens = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/pathogen_list.txt")
        filtered_data = process_annotation(data_to_be_processed, mge_len_file, pathogens)
        result = calculate_score(myargs['-c'], filtered_data, pipeline)
	elif pipeline == 2:        
        data_to_be_processed = [annotated_data[3], annotated_data[1], annotated_data[2]]
        pathogens = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/eskape.txt")
        filtered_data = process_annotation(data_to_be_processed, mge_len_file, pathogens)
        result = calculate_score(myargs['-c'], filtered_data, pipeline)
    else:
        # Eco Risk
        data_to_be_processed_e = [annotated_data[0], annotated_data[1], annotated_data[2]]
        pathogens_e = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/pathogen_list.txt")
        filtered_data_e = process_annotation(data_to_be_processed_e, mge_len_file, pathogens_e)
        result_e = calculate_score(myargs['-c'], filtered_data_e, 1)
        
        # Human Risk
        data_to_be_processed_h = [annotated_data[3], annotated_data[1], annotated_data[2]]
        pathogens_h = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/eskape.txt")
        filtered_data_h = process_annotation(data_to_be_processed_h, mge_len_file, pathogens_h)
        result_h = calculate_score(myargs['-c'], filtered_data_h, 2)
        
        result = pd.concat([result_e, result_h])
    
    result.to_csv(out_file, header=True, index=None, sep="\t")
    print(f"Results saved to {out_file}")