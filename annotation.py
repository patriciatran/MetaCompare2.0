import os, sys, glob
import subprocess
import pandas as pd
import shutil

import warnings
warnings.simplefilter(action='ignore', category=Warning)

IDENTITY = 60
MIN_ALIGN_LENGTH = 25
COVG_OF_ALIGN_MGE = 0.9

def generate_annotation(contig_file, out_dir, nthread = '64'):
	sample_name = contig_file.split('/')[-1].split('.')[0]
	prodigal_file = os.path.join(out_dir, sample_name + "_genes." + contig_file.split('/')[-1].split('.')[1]) 
	prodigal_file_faa = os.path.join(out_dir, sample_name + "_proteins." + contig_file.split('/')[-1].split('.')[1]) 

	if not os.path.exists(prodigal_file):
		print("Running Prodigal")
		subprocess.call(["pprodigal", "-i", contig_file, "-d", prodigal_file, "-a", prodigal_file_faa, "-p", "meta", "-C", "4000", "-T", "32", "-o", os.path.join(out_dir, "prodigal_log")])
	else:
		print("Skipping: Prodigal output already exists")

	
	arg_name = sample_name + "_ARG.csv"		
	if not os.path.exists(os.path.join(out_dir, arg_name)):
		print('Running Diamond Blastx on ARGDB')
		subprocess.call(["diamond", "blastx", "-d", os.path.dirname(os.path.abspath(__file__))+"/metacmpDB/ARGDB", "--query", prodigal_file, \
						"--out", os.path.join(out_dir, arg_name), "--outfmt", "6", \
						"--threads", nthread, "--evalue", "1e-10"])
	else:
		print('Skipping: Diamond output against ARGs already exists')
		
	arg_name_hh = sample_name + "_hh_ARG.csv"	
	if not os.path.exists(os.path.join(out_dir, arg_name_hh)):
		print('Running Diamond Blastx on ARGDB_hh')
		subprocess.call(["diamond", "blastx", "-d", os.path.dirname(os.path.abspath(__file__))+"/metacmpDB/ARGDB_hh", "--query", prodigal_file, \
						"--out", os.path.join(out_dir, arg_name_hh), "--outfmt", "6", \
						"--threads", nthread, "--evalue", "1e-10"])
	else:
		print('Skipping: Diamond output against ARGs already exists')
		
	mge_name = sample_name + "_MGE.csv"
	if not os.path.exists(os.path.join(out_dir, mge_name)):
		print('Running Diamond Blastx on MGEDB')
		subprocess.call(["diamond", "blastx", "-d", os.path.dirname(os.path.abspath(__file__))+"/metacmpDB/MGEDB", "--query", prodigal_file, \
						"--out", os.path.join(out_dir, mge_name), "--outfmt", "6", \
						"--threads", nthread, "--evalue", "1e-10"])
	else:
		print('Skipping: Diamond output against MGEs already exists')

	pathogen_name = sample_name + "_Pathogens.tsv"
	if not os.path.exists(os.path.join(out_dir, pathogen_name)):
		print('Running mmseq2 on GTDB')
		command = ['sh', os.path.dirname(os.path.abspath(__file__))+"mmseq.sh", contig_file, out_dir, pathogen_name]		
		if(out_dir != "" and not out_dir.endswith('/')):
			out_dir += "/"		
		subprocess.call(['sh', os.path.dirname(os.path.abspath(__file__))+"/./mmseq.sh", contig_file, out_dir, pathogen_name])
		try:
			shutil.rmtree(os.path.join(out_dir, "sample.tmpFolder"))
			for f in glob.glob(os.path.join(out_dir, "sample.contigs*")):
				os.remove(f)
			for f in glob.glob(os.path.join(out_dir, "sample.assignments*")):
				os.remove(f)
		except OSError as e:
			print("Error: %s - %s." % (e.filename, e.strerror))
	else:
		print('Skipping: mmseq2 output against Pathogens already exists')
	
	return [os.path.join(out_dir, arg_name), os.path.join(out_dir, mge_name), os.path.join(out_dir, pathogen_name), os.path.join(out_dir, arg_name_hh)]

def filter_diamond(filename):
	data = pd.read_csv(filename, sep='\t', header=None)
	data.columns = ['id', 'sub_id', 'identity', 'alignLen', 'mismat', 'gapOpens', 'qStart', 'qEnd', 'sStart', 'sEnd', 'eval', 'bit']
		
	# filter out contigs identity under 60
	iden_filtered = data[data.identity > IDENTITY]
		
	# filter out contigs alignment length under 25
	filtered_data = iden_filtered[iden_filtered.alignLen > MIN_ALIGN_LENGTH]
	
	return filtered_data

def process_annotation(data_files, mge_len_file, pathogen_list):	
	arg_file = data_files[0]
	mge_file = data_files[1]
	path_file = data_files[2]
	#arg_hh_file = data_files[3]
	
	# Open blast output against ARG DB
	if not os.path.getsize(arg_file) > 0:
		#file is empty
		print('Warning: '+ arg_file+ ' is empty.')
		arg_data = pd.DataFrame()
	else:		
		arg_data = filter_diamond(arg_file)		
		# Note: 'arg_data' data frame contains name and position of Antibiotic Resistence Genes in contigs		
		
	# Open blast output against MGE DB
	if not os.path.getsize(mge_file) > 0:
		#file is empty
		print('Warning: '+ mge_file + ' is empty.')
		mge_data = pd.DataFrame()
	else:
		mge_data = filter_diamond(mge_file) 
		# filter out contigs having less than 90% coverage of the reference
		if not os.path.exists(mge_len_file):
			print(mge_len_file + ' file does not exists.')
			sys.exit()
		else:
			mge_len = pd.read_csv(mge_len_file, sep='\t', header=None)
			mge_len.columns = ['sub_id', 'ref_gene_leng']
			mge_merged = pd.merge(mge_data, mge_len, how = 'left', on = 'sub_id')
			mge_final = mge_merged[mge_merged.alignLen > (mge_merged.ref_gene_leng * COVG_OF_ALIGN_MGE)]		
			# Note: 'mge_final' data frame contains name and position of MGEs in contigs
            
	# Open mmseq output against GTDB database 
	if not os.path.getsize(path_file) > 0:
		#file is empty 
		print('Warning: '+ path_file + ' is empty.')
		path_all = pd.DataFrame()
	else:
		path_data = pd.read_csv(path_file, sep='\t', header=None)
		path_data.columns = ['id', 'NCBI_ID','rank', 'name', 'nPass', 'nRetain', 'nAssign', 'bit', 'taxonomy']
		path_filtered = path_data.loc[path_data['rank'].isin(["family", "genus", "species", "strain"])]
		
		def get_pathogens(path_data, path_file):
			pathogens = pd.read_csv(path_file, sep = "\t")
			ranks = pathogens['rank'].drop_duplicates().values.tolist()
			path_final = pd.DataFrame(columns = path_data.columns)
			for rank in ranks:
				filtered_data = path_data.loc[path_data['rank'].isin([rank])]
				path_of_rank = pathogens.loc[pathogens['rank'] == rank]["name"].values.tolist()
				selected_data=filtered_data.loc[filtered_data['name'].isin(path_of_rank)]
				path_final = pd.concat([path_final, selected_data])
			
			### only for strain
			filtered_data = path_data.loc[path_data['rank'].isin(["strain"])]
			if not filtered_data.empty:
				filtered_data['name'] = filtered_data['name'].str.split(' ').str[:2]
				filtered_data['name'] = filtered_data['name'].apply(lambda x: ' '.join(map(str, x)))
				path_of_rank = pathogens.loc[pathogens['rank'] == "species"]["name"].values.tolist()
				selected_data=filtered_data.loc[filtered_data['name'].isin(path_of_rank)]
				path_final = pd.concat([path_final, selected_data])
			return path_final
				
		path_all = get_pathogens(path_filtered, pathogen_list)
		
	return [arg_data, mge_final, path_all]
