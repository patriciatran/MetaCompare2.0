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
        print('\nUsage: ./metacompare.py -c filename1.fa [-t 64 -b 1 -s path_prefix -db path] \n')
        print('\t-c: Specify FASTA file containing assembled contigs [required]')
        print('\t-t: Specify the number of threads (default: 64).')
        print('\t-b: Specify pipeline [0: both (default), 1: eco risk, 2: human risk].')
        print('\t-o: Output file path.')
        print('\t-s: Skip annotation. Provide the prefix/path to existing Prodigal files.')
        print('\t-db: Path to the GTDB database for mmseqs.')
        print()
        exit()

    # If the database path is not provided for the mmseq.sh script, we point to the old default location (which should be the same folder location as the data)
    default_db = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/GTDB/gtdb")
    db_path = myargs.get('-db', default_db)

    nthread = myargs.get('-t', '64')
    out_dir = myargs.get('-o', '')
    pipeline = int(myargs.get('-b', '0'))
    
    mge_len_file = os.path.join(os.path.dirname(os.path.abspath(__file__)), "metacmpDB/MGE_len.txt")
    sample_name = os.path.splitext(os.path.basename(myargs['-c']))[0]
    out_file = os.path.join(out_dir, sample_name + "_out.txt")    
    
    if '-s' in myargs:
        print(f"Skipping annotation. Using prefix: {myargs['-s']}")
        prefix = myargs['-s']
        annotated_data = [
            f"{prefix}.gene.faa", 
            f"{prefix}.gene.fna", 
            f"{prefix}.gff",      
            f"{prefix}.renamed"   
        ]
    else:
        print(f"Running annotation step using database at: {db_path}")
        # 2. Pass db_path to the annotation function
        # NOTE: You must update the signature of generate_annotation in annotation.py to accept this!
        annotated_data = generate_annotation(myargs['-c'], out_dir, nthread, db_path)    
    
    # Pipeline execution...
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