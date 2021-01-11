'''
This program allows you to score the robustness (or lack-thereof) of a protein sequence
for which you have estimates of the functional effects of most amino acid substitutions.

Run with the following call:

python Robustness_scores.py --setup [setup_file.csv]

The setup file is organized into a Notes section, which is ignored, and a Jobs section
that lets you program how you want the script to run. It is formatted as follows
(shown tab-delimited, but it is .csv for use with Excel)

Scores_file						DNA_sequence_file	seq_ID	Species	Analysis_type	Iterations	Filter_sites	Codon_bias	Transform_matrix	output_file
Standard_protein_set/Brca1.tsv	Robu_genes.fasta	Brca1	human	Sequence		100			TRUE			No bias							output/Brca1_test.tsv
Standard_protein_set/Erk.tsv	Robu_genes.fasta	Brca1	human	Position					False			No bias							

The header should not be altered. However, columns that are not required can be left empty.
Scores_file, DNA_sequence_file, seq_ID, and Analysis_type are required with no default.
Descriptions of the fields are in the Notes section.

The script takes as input, a tab-delimited Scores file with the following format:

position	AAwt	AAmut	norm_score
2	A	C	0.814995277123
2	A	D	0.926545681819
2	A	E	1.07835366716
2	A	F	0.830394119876
...

These exact column names are required, and additional columns will be ignored. Scores
should be normalized such that the wild type protein has a norm_score of 1 and 'dead',
'null', or premature termination variants have a norm_score of 0. Scores >1 or <0 will
be considered according the values specified, but may not represent a real functional
effect, and extreme outliers can drive a strong effect that may have more to do with
technical noise than with real biology. In the paper, we compress all scores >2 or <-1
to 2 or -1, respectively, although you may consider using different boundaries.

As an additional input, you must provide a DNA fasta file for the protein sequence 
under consideration. The sequence in the fasta should start and stop where the scores
file starts and stops. Gaps in the scores file can be filled in in the fasta with the
real DNA sequence or with dashes; they should not be clipped out. The position column
in the scores file does not need to correspond with the position in the fasta sequence.
The script will attempt to line up the sequence in the fasta to the sequence it infers
from your scores file. If it does not match, it will inform you and exit.

This script should process 1000 iterations in Sequence mode in around a minute on a 
desktop for a typical-length protein, and run in Position mode in a few seconds. For our
longest dataset, 757 positions, the script ran 10,000 iterations in 17 minutes. It
will exit with a warning if it detects certain invalid inputs. Therefore, if you plan
to run many iterations (ie. >10,000), it is recommended to run it first with all of the
settings specified, but set iterations to 10 to check that no errors are thrown and it
processes all of the jobs. Then the iterations can be increased to any reasonable number.

The script runs all iterations for a single job before it starts printing, so the number
of iterations that can be run in one call is memory-limited. We have run 100,000
iterations without memory issues, but going much above that, it would be recommended to
break it into multiple jobs with different output files, then concatenate the output
files.
'''

def main(setup_file):
	
	jobs_to_do = set_up_jobs(setup_file)
	print(jobs_to_do)
	for i in range(1,len(jobs_to_do)+1):
		run_job(jobs_to_do[i])
	
def run_job(params):		
	#Initialize various objects
	
	scores_tsv = params[0]
	DNA_file = params[1]
	protein_name = params[2]
	species = params[3]
	analysis = params[4]
	iterations = params[5]
	filter = params[6]
	codon_bias = params[7]
	transformation_matrix_file = params[8]
	output_file = params[9]
	
	if len(protein_name) == 0:
		protein_name = scores_tsv.split(".")[-2].split("/")[-1]
	wt_DNA_seq = read_fasta(DNA_file)[protein_name]
	wt_pro_seq = translate(wt_DNA_seq)
	start_pos, temp_pro_seq, end_pos = infer_protein_seq(scores_tsv)
	scores_dict = read_scores_tsv(scores_tsv) #keys are tuples with structure (int(position), str(AA_wt), str(AA_mut)), values are 0-1 scaled log enrichments

	positions_to_evaluate = []
	output_dict = {}
	
	if len(transformation_matrix_file) > 0:
		transform_mat = read_transformation_matrix(file_name = transformation_matrix_file)
	else:
		transform_mat = None
	
	#Filtering and checks
	#Compare protein sequence from translating the DNA to the protein sequence inferable from scores data
	if match_protein_sequences(temp_pro_seq, wt_pro_seq) == False:
		sys.stderr.write("Error: DNA sequence provided does not match that inferred from scores_tsv\nFrom Fasta:\t" + temp_pro_seq + '\nFrom Scores:\t' + wt_pro_seq + '\n')
		return None
	#Impose filtering on sites (if specified) and throw out sites with missing data for any neighbor of a synonymous codon
	if filter.upper() == "TRUE":
		for i in range(start_pos, end_pos+1):
			if check_if_position_complete(i, scores_dict, start_pos, wt_DNA_seq) == True:
				positions_to_evaluate.append(i)
	else: #As of 1/2/19, this no-filter version is supposed to be functional. CRS_score() will probably throw a divide-by-zero error if no scores at all are available for position, but this hasn't been tested
		for i in range(start_pos, end_pos+1):
			positions_to_evaluate.append(i)
	sys.stderr.write("Data loaded successfully and analysis started for protein: "+protein_name+"\n")
	
	#Analysis time!
	
	#First, go through the protein and define the average robustness for all synonymous codons at each position
	
	Avg_CRS_dict = {}
	for i in positions_to_evaluate:
		CRS_dict = {}
		wt_codon = wt_DNA_seq[3*(i - start_pos):3*(i - start_pos)+3]
		wt_aa = translation_table[wt_codon]
		CRS_dict[wt_codon] = CRS_score(wt_DNA_seq, i, scores_dict, start_pos, tra_mat = transform_mat, soo = species)
		for c in codon_synonyms[wt_codon]:
			CRS_dict[c] = CRS_score(wt_DNA_seq, i, scores_dict, start_pos, set_codon = c, tra_mat = transform_mat, soo = species)
		Avg_CRS_dict[i] = sum(CRS_dict.values()) / float(len(CRS_dict))
	
	#Then, analysis varies by mode
	#For sequence, go through codons and add CRS to wt_CRS
	#Then iteratively do the same for syn_CRS
	if analysis == "Sequence":
		iterations = int(iterations)
		#Standard analysis for comparing the full-length wt sequence against many full-length synonymous sequences
		wt_CRS = []
		for i in positions_to_evaluate:
			wt_CRS.append(CRS_score(wt_DNA_seq, i, scores_dict, start_pos, tra_mat = transform_mat, soo = species) - Avg_CRS_dict[i])
		SRS_wt = sum(wt_CRS) / float(len(wt_CRS))
		
		#Iteratively generate synonymous sequences, score them, and compare them to wt -- write values to output_dict
		for j in range(0,iterations):
			if codon_bias.upper() == "TRUE":
				syn_DNA_seq = reverse_translate(wt_pro_seq, species)
			
			else:
				syn_DNA_seq = reverse_translate(wt_pro_seq)
			syn_CRS = []
			for i in positions_to_evaluate:
				syn_CRS.append(CRS_score(syn_DNA_seq, i, scores_dict, start_pos, tra_mat = transform_mat, soo = species) - Avg_CRS_dict[i])
			SRS_syn = sum(syn_CRS) / float(len(syn_CRS))
			output_dict[j] = '\t'.join([str(SRS_wt), str(SRS_syn), str(SRS_wt / SRS_syn), str(SRS_wt - SRS_syn)])
			if j > 0 and j % 200 == 0:
				sys.stderr.write("Completed %d iterations\n" % j) #Importantly, all analysis is done before printing, this should help assess runtime
		#Print things
		if len(output_file) > 4:
			file_path = output_file
			print(file_path)
		else:		
			file_path = "./output/" + '_'.join([protein_name,analysis,iterations]) + ".tsv"
		directory = os.path.dirname(file_path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		with open(file_path,'w+') as file_out:
			file_out.write('\t'.join(['avg_wt', 'avg_syn', 'avg_ratio', 'avg_diff']))
			for j in range(0,iterations):
				file_out.write('\n'+output_dict[j])
		sys.stderr.write("Completed without error, written to " + file_path + '\n')
	
	elif analysis == "Position":
		#For each codon in positions_to_evaluate, calculate CRS for wt and all synonyms, then take average over all codons at position. For each codon, print CRS and comparison to average (report whether its wt or syn too)
		output_dict = {0:["Protein","Position","AAwt","codon","wt?","CRS","CRS_avg","diff"]}
		for i in positions_to_evaluate:
			CRS_dict = {}
			wt_codon = wt_DNA_seq[3*(i - start_pos):3*(i - start_pos)+3]
			wt_aa = translation_table[wt_codon]
			CRS_dict[wt_codon] = CRS_score(wt_DNA_seq, i, scores_dict, start_pos, tra_mat = transform_mat, soo = species)
			for c in codon_synonyms[wt_codon]:
				CRS_dict[c] = CRS_score(wt_DNA_seq, i, scores_dict, start_pos, set_codon = c, tra_mat = transform_mat, soo = species)
			CRS_avg = sum(CRS_dict.values()) / float(len(CRS_dict))
			for c in CRS_dict:
				if c == wt_codon:
					w = "wt"
				else:
					w = "syn"
				output_dict[len(output_dict)] = [protein_name,str(i),wt_aa,c,w,str(CRS_dict[c]),str(CRS_avg),str(CRS_dict[c] - CRS_avg)]
		
		#Print things
		file_path = "./output/" + protein_name + ".out1"
		directory = os.path.dirname(file_path)
		if not os.path.exists(directory):
			os.makedirs(directory)
		with open(file_path,'w+') as file_out:
			for j in range(0,len(output_dict)):
				file_out.write('\t'.join(output_dict[j]) + '\n')
		sys.stderr.write("Completed without error, written to " + file_path + '\n')
	
	elif analysis == 2:
		pass
		'''
            Fill this in later
            calculate CRS for every codon in iterated syn seqs
            This should be the 'rawest' output
		'''	
	else:
		sys.stderr.write("problem encountered. Probably the analysis type specified was not valid. Use --analysis 0/1/2, or leave off flag for analysis 0")

def set_up_jobs(setup_file):
	with open(setup_file,'r') as csv:
		jobs = {}
		read_jobs = False
		job_count = 0
		column_indices = {}
		for line in csv:
			if read_jobs == True:
				j = line.strip().split(',')
				if len(j) > 7:
					job_count +=1
				if job_count == 1:
					for i in range(0,len(j)):
						column_indices[j[i]] = i
				elif job_count > 1:
					jobs[job_count - 1] = (j[column_indices["Scores_file"]], j[column_indices["DNA_sequence_file"]], j[column_indices["seq_ID"]], j[column_indices["Species"]], j[column_indices["Analysis_type"]], j[column_indices["Iterations"]], j[column_indices["Filter_sites"]], j[column_indices["Codon_bias"]], j[column_indices["Transform_matrix"]], j[column_indices["output_file"]])
					
			elif line.startswith("#Jobs"):
				read_jobs = True
	return jobs	
			
def check_if_position_complete(position, scores_dict, start_pos, wt_DNA_seq):
	check = True
	wt_codon = wt_DNA_seq[(position - start_pos)*3:(position - start_pos)*3+3]
	for j in codon_neighbors[wt_codon]:
		if translation_table[wt_codon] != translation_table[j] and translation_table[j] != '*':
			sub_name = (int(position), translation_table[wt_codon], translation_table[j])
			#print sub_name, scores_dict.get(sub_name)
			if scores_dict.get(sub_name) == None:
				check = False			
	for syn_codon in codon_synonyms[wt_codon]:
		for k in codon_neighbors[syn_codon]:
			if translation_table[syn_codon] != translation_table[k] and translation_table[k] != '*' and translation_table[syn_codon] != '*':
				sub_name = (int(position), translation_table[syn_codon], translation_table[k])
				if scores_dict.get(sub_name) == None:
					check = False
	if check == False:
		sys.stderr.write("position %s excluded due to missing data\n" % position)
	return check
	
def infer_protein_seq(scores):
	with open(scores,'r') as tsv:
		line_count = 0
		seq_dict = {}
		for line in tsv:
			#print line
			row = line.strip().split('\t')
			#print '\t'.join(row)
			line_count += 1
			if line_count == 1:
				wt_ind = row.index("AAwt")
				pos_ind = row.index("position")
			elif seq_dict.get(row[pos_ind]) == None:
				seq_dict[int(row[pos_ind])] = row[wt_ind]
		#print seq_dict
	tsv.close()
	seq = []
	for k in range(min(seq_dict.keys()),max(seq_dict.keys())+1):
		seq.append(seq_dict.get(k,"-"))
	return min(seq_dict.keys()), ''.join(seq), max(seq_dict.keys()) 

def read_fasta(file):
	seqs = {}
	with open(file,'r') as fasta:
		for line in fasta:
			if line.startswith(">"):
				header = line.strip()[1:]
				seqs[header] = ""
			else:
				seqs[header] += line.strip()
	return seqs

def translate(DNA_sequence, frame = 0):
	DNA_sequence = DNA_sequence.upper()
	protein_sequence_list = []
	for i in range(frame,len(DNA_sequence),3):
		protein_sequence_list.append(translation_table[DNA_sequence[i:i+3]])
	return ''.join(protein_sequence_list)

def reverse_translate(protein, species = "neutral"):
	DNA = []
	for aa in protein:
		if species == "neutral":
			DNA.append(numpy.random.choice(codon_freqs[aa]['codons']).upper())
		else:
			DNA.append(numpy.random.choice(codon_freqs[aa]['codons'], p = codon_freqs[aa][species]).upper())
	return ''.join(DNA)

def match_protein_sequences(seq1, seq2):
	if len(seq1) != len(seq2):
		sys.stderr.write("Error: Protein sequences don't match length between scores file and fasta\nFrom Fasta:\t")
		sys.stderr.write(seq2+'\nFrom Scores:\t')
		sys.stderr.write(seq1+'\n')
		sys.exit()
	else:
		test = True
		for i in range(0,len(seq1)):
			if seq1[i] != '-' and seq2[i] != '-' and seq1[i] != seq2[i]:
				test = False
	return test

def read_scores_tsv(file_name):
	All_subs = {}
	line_count = 0
	with open(file_name,'r') as scores:
		for line in scores:
			line_count += 1
			row = line.strip().split('\t')
			if line_count == 1:
				features = row
				pos_ind = features.index("position")
				aa1_ind = features.index("AAwt")
				aa2_ind = features.index("AAmut")
				score_ind = features.index("norm_score")
			else:
				sub_name = (int(row[pos_ind]), row[aa1_ind], row[aa2_ind])
				All_subs[sub_name] = float(row[score_ind])
	return All_subs

def read_transformation_matrix(file_name = None, species = None, process = None):
	if file_name == None:
		if species == None or process == None:
			sys.stderr.write("failed to find transformation matrix\nexiting\n")
			sys.exit()
		else:
			file_name = "./" + species.lower() + "_" + process.lower() + "_transform.txt"
		
	if not os.path.exists(file_name):
		sys.stderr.write("failed to find transformation matrix\nexiting\n")
		sys.exit()
	transform_matrix = {}
	with open(file_name,'r') as mat_file:
		Block = None
		features = {}
		for line in mat_file:
			if line.startswith(">"):
				Block = line.strip()[1:]
			elif Block == "General_info":
				if line.count(':') == 1:
					features[line.strip().split(':')[0]] = line.strip().split(':')[1]
			
			elif Block == "Array":
				if line.strip().split('\t')[0] == "codon":
					cols = {}
					for i in range(1,65):
						cols[i] = line.strip().split('\t')[i]
				elif line[0] in ['A','C','T','G']:
					row = line.strip().split('\t')
					transform_matrix[row[0]] = {}
					for i in range(1,65):
						transform_matrix[row[0]][cols[i]] = float(row[i])
		mat_file.close()
	return transform_matrix			
	
	
def CRS_score(DNA_seq, position, scores_table, start_pos, tra_mat = None, set_codon = None, soo = "neutral"): # soo being species-of-origin, with name 'species' being used only as a main() variable to avoid confusion, but they mean the same thing
	if tra_mat == None:
		tra_mat = neutral_transformation_matrix
	elif soo == "neutral":
		sys.stderr.write("Error: if a weight matrix is specified, the species must also be specified")
		sys.exit()
	i = 3*(position - start_pos)
	if set_codon == None:
		codon = DNA_seq[i:i+3]
	else:
		codon = set_codon
	sum = 0.0
	denom = 0.0
	for sub in sorted(codon_neighbors.keys()):
		if translation_table[codon] == translation_table[sub]:
			score = 1.0
		elif translation_table[sub] == '*':
			score = 0.0
		else:
			mut = (position, translation_table[codon], translation_table[sub])
			try:
				score = float(scores_table[mut])
			except KeyError:
				score = "missing"
		if tra_mat == None:
			weight = 1
		else:
			weight = tra_mat[codon][sub]	
		if score != "missing":
			sum += (score * weight)
			denom += weight
	CRS = sum / denom
	return CRS

translation_table = {	
	'ACC': "T", 'ATG': "M", 'ACA': "T", 
	'ACG': "T", 'ATC': "I", 'AAC': "N", 
	'ATA': "I", 'AGG': "R", 'CCT': "P", 
	'CTC': "L", 'AGC': "S", 'AAG': "K", 
	'AGA': "R", 'CAT': "H", 'AAT': "N", 
	'ATT': "I", 'CTG': "L", 'CTA': "L", 
	'ACT': "T", 'CAC': "H", 'AAA': "K", 
	'CCG': "P", 'AGT': "S", 'CCA': "P", 
	'CAA': "Q", 'CCC': "P", 'TAT': "Y", 
	'GGT': "G", 'TGT': "C", 'CGA': "R", 
	'CAG': "Q", 'CGC': "R", 'GAT': "D", 
	'CGG': "R", 'CTT': "L", 'TGC': "C", 
	'GGG': "G", 'TAG': "*", 'GGA': "G", 
	'TAA': "*", 'GGC': "G", 'TAC': "Y", 
	'GAG': "E", 'TCG': "S", 'TTT': "F", 
	'GAC': "D", 'CGT': "R", 'GAA': "E", 
	'TCA': "S", 'GCA': "A", 'GTA': "V", 
	'GCC': "A", 'GTC': "V", 'GCG': "A", 
	'GTG': "V", 'TTC': "F", 'GTT': "V", 
	'GCT': "A", 'TTA': "L", 'TGA': "*", 
	'TTG': "L", 'TCC': "S", 'TGG': "W", 
	'TCT': "S"}
	
codon_neighbors = {	
	'ACC': ['CCC', 'GCC', 'TCC', 'AAC', 'AGC', 'ATC', 'ACA', 'ACG', 'ACT'], 
	'ATG': ['CTG', 'GTG', 'TTG', 'ACG', 'AGG', 'AAG', 'ATC', 'ATA', 'ATT'], 
	'ACA': ['CCA', 'GCA', 'TCA', 'AAA', 'AGA', 'ATA', 'ACC', 'ACG', 'ACT'], 
	'ACG': ['CCG', 'GCG', 'TCG', 'AAG', 'AGG', 'ATG', 'ACC', 'ACA', 'ACT'], 
	'ATC': ['CTC', 'GTC', 'TTC', 'ACC', 'AGC', 'AAC', 'ATA', 'ATG', 'ATT'], 
	'AAC': ['CAC', 'GAC', 'TAC', 'ACC', 'AGC', 'ATC', 'AAA', 'AAG', 'AAT'], 
	'ATA': ['CTA', 'GTA', 'TTA', 'ACA', 'AGA', 'AAA', 'ATC', 'ATG', 'ATT'], 
	'AGG': ['CGG', 'GGG', 'TGG', 'ACG', 'AAG', 'ATG', 'AGC', 'AGA', 'AGT'], 
	'CCT': ['ACT', 'GCT', 'TCT', 'CAT', 'CGT', 'CTT', 'CCC', 'CCG', 'CCA'], 
	'ACT': ['CCT', 'GCT', 'TCT', 'AAT', 'AGT', 'ATT', 'ACC', 'ACG', 'ACA'], 
	'AGC': ['CGC', 'GGC', 'TGC', 'ACC', 'AAC', 'ATC', 'AGA', 'AGG', 'AGT'], 
	'AAG': ['CAG', 'GAG', 'TAG', 'ACG', 'AGG', 'ATG', 'AAC', 'AAA', 'AAT'], 
	'AGA': ['CGA', 'GGA', 'TGA', 'ACA', 'AAA', 'ATA', 'AGC', 'AGG', 'AGT'], 
	'CAT': ['AAT', 'GAT', 'TAT', 'CCT', 'CGT', 'CTT', 'CAC', 'CAG', 'CAA'], 
	'AAT': ['CAT', 'GAT', 'TAT', 'ACT', 'AGT', 'ATT', 'AAC', 'AAG', 'AAA'], 
	'ATT': ['CTT', 'GTT', 'TTT', 'ACT', 'AGT', 'AAT', 'ATC', 'ATG', 'ATA'], 
	'CTG': ['ATG', 'GTG', 'TTG', 'CCG', 'CGG', 'CAG', 'CTC', 'CTA', 'CTT'], 
	'CTA': ['ATA', 'GTA', 'TTA', 'CCA', 'CGA', 'CAA', 'CTC', 'CTG', 'CTT'], 
	'CTC': ['ATC', 'GTC', 'TTC', 'CCC', 'CGC', 'CAC', 'CTA', 'CTG', 'CTT'], 
	'CAC': ['AAC', 'GAC', 'TAC', 'CCC', 'CGC', 'CTC', 'CAA', 'CAG', 'CAT'], 
	'AAA': ['CAA', 'GAA', 'TAA', 'ACA', 'AGA', 'ATA', 'AAC', 'AAG', 'AAT'], 
	'CCG': ['ACG', 'GCG', 'TCG', 'CAG', 'CGG', 'CTG', 'CCC', 'CCA', 'CCT'], 
	'AGT': ['CGT', 'GGT', 'TGT', 'ACT', 'AAT', 'ATT', 'AGC', 'AGG', 'AGA'], 
	'CCA': ['ACA', 'GCA', 'TCA', 'CAA', 'CGA', 'CTA', 'CCC', 'CCG', 'CCT'], 
	'CAA': ['AAA', 'GAA', 'TAA', 'CCA', 'CGA', 'CTA', 'CAC', 'CAG', 'CAT'], 
	'CCC': ['ACC', 'GCC', 'TCC', 'CAC', 'CGC', 'CTC', 'CCA', 'CCG', 'CCT'], 
	'TAT': ['CAT', 'GAT', 'AAT', 'TCT', 'TGT', 'TTT', 'TAC', 'TAG', 'TAA'], 
	'GGT': ['CGT', 'AGT', 'TGT', 'GCT', 'GAT', 'GTT', 'GGC', 'GGG', 'GGA'], 
	'TGT': ['CGT', 'GGT', 'AGT', 'TCT', 'TAT', 'TTT', 'TGC', 'TGG', 'TGA'], 
	'CGA': ['AGA', 'GGA', 'TGA', 'CCA', 'CAA', 'CTA', 'CGC', 'CGG', 'CGT'], 
	'CAG': ['AAG', 'GAG', 'TAG', 'CCG', 'CGG', 'CTG', 'CAC', 'CAA', 'CAT'], 
	'CGC': ['AGC', 'GGC', 'TGC', 'CCC', 'CAC', 'CTC', 'CGA', 'CGG', 'CGT'], 
	'GAT': ['CAT', 'AAT', 'TAT', 'GCT', 'GGT', 'GTT', 'GAC', 'GAG', 'GAA'], 
	'CGG': ['AGG', 'GGG', 'TGG', 'CCG', 'CAG', 'CTG', 'CGC', 'CGA', 'CGT'], 
	'CTT': ['ATT', 'GTT', 'TTT', 'CCT', 'CGT', 'CAT', 'CTC', 'CTG', 'CTA'], 
	'TGC': ['CGC', 'GGC', 'AGC', 'TCC', 'TAC', 'TTC', 'TGA', 'TGG', 'TGT'], 
	'GGG': ['CGG', 'AGG', 'TGG', 'GCG', 'GAG', 'GTG', 'GGC', 'GGA', 'GGT'], 
	'TAG': ['CAG', 'GAG', 'AAG', 'TCG', 'TGG', 'TTG', 'TAC', 'TAA', 'TAT'], 
	'GGA': ['CGA', 'AGA', 'TGA', 'GCA', 'GAA', 'GTA', 'GGC', 'GGG', 'GGT'], 
	'TAA': ['CAA', 'GAA', 'AAA', 'TCA', 'TGA', 'TTA', 'TAC', 'TAG', 'TAT'], 
	'GGC': ['CGC', 'AGC', 'TGC', 'GCC', 'GAC', 'GTC', 'GGA', 'GGG', 'GGT'], 
	'TAC': ['CAC', 'GAC', 'AAC', 'TCC', 'TGC', 'TTC', 'TAA', 'TAG', 'TAT'], 
	'GAG': ['CAG', 'AAG', 'TAG', 'GCG', 'GGG', 'GTG', 'GAC', 'GAA', 'GAT'], 
	'TCG': ['CCG', 'GCG', 'ACG', 'TAG', 'TGG', 'TTG', 'TCC', 'TCA', 'TCT'], 
	'TTT': ['CTT', 'GTT', 'ATT', 'TCT', 'TGT', 'TAT', 'TTC', 'TTG', 'TTA'], 
	'GAC': ['CAC', 'AAC', 'TAC', 'GCC', 'GGC', 'GTC', 'GAA', 'GAG', 'GAT'], 
	'CGT': ['AGT', 'GGT', 'TGT', 'CCT', 'CAT', 'CTT', 'CGC', 'CGG', 'CGA'], 
	'GAA': ['CAA', 'AAA', 'TAA', 'GCA', 'GGA', 'GTA', 'GAC', 'GAG', 'GAT'], 
	'TCA': ['CCA', 'GCA', 'ACA', 'TAA', 'TGA', 'TTA', 'TCC', 'TCG', 'TCT'], 
	'GCA': ['CCA', 'ACA', 'TCA', 'GAA', 'GGA', 'GTA', 'GCC', 'GCG', 'GCT'], 
	'GTA': ['CTA', 'ATA', 'TTA', 'GCA', 'GGA', 'GAA', 'GTC', 'GTG', 'GTT'], 
	'GCC': ['CCC', 'ACC', 'TCC', 'GAC', 'GGC', 'GTC', 'GCA', 'GCG', 'GCT'], 
	'GTC': ['CTC', 'ATC', 'TTC', 'GCC', 'GGC', 'GAC', 'GTA', 'GTG', 'GTT'], 
	'GCG': ['CCG', 'ACG', 'TCG', 'GAG', 'GGG', 'GTG', 'GCC', 'GCA', 'GCT'], 
	'GTG': ['CTG', 'ATG', 'TTG', 'GCG', 'GGG', 'GAG', 'GTC', 'GTA', 'GTT'], 
	'TTC': ['CTC', 'GTC', 'ATC', 'TCC', 'TGC', 'TAC', 'TTA', 'TTG', 'TTT'], 
	'GTT': ['CTT', 'ATT', 'TTT', 'GCT', 'GGT', 'GAT', 'GTC', 'GTG', 'GTA'], 
	'GCT': ['CCT', 'ACT', 'TCT', 'GAT', 'GGT', 'GTT', 'GCC', 'GCG', 'GCA'], 
	'TTA': ['CTA', 'GTA', 'ATA', 'TCA', 'TGA', 'TAA', 'TTC', 'TTG', 'TTT'], 
	'TGA': ['CGA', 'GGA', 'AGA', 'TCA', 'TAA', 'TTA', 'TGC', 'TGG', 'TGT'], 
	'TTG': ['CTG', 'GTG', 'ATG', 'TCG', 'TGG', 'TAG', 'TTC', 'TTA', 'TTT'], 
	'TCC': ['CCC', 'GCC', 'ACC', 'TAC', 'TGC', 'TTC', 'TCA', 'TCG', 'TCT'], 
	'TGG': ['CGG', 'GGG', 'AGG', 'TCG', 'TAG', 'TTG', 'TGC', 'TGA', 'TGT'], 
	'TCT': ['CCT', 'GCT', 'ACT', 'TAT', 'TGT', 'TTT', 'TCC', 'TCG', 'TCA']}

codon_synonyms = {
	'CTT': ['CTG', 'CTA', 'CTC', 'TTA', 'TTG'],
	'ATG': [],
	'AAG': ['AAA'],
	'AAA': ['AAG'],
	'ATC': ['ATA', 'ATT'],
	'AAC': ['AAT'],
	'ATA': ['ATC', 'ATT'],
	'AGG': ['AGA', 'CGA', 'CGC', 'CGG', 'CGT'],
	'CCT': ['CCG', 'CCA', 'CCC'],
	'CTC': ['CTG', 'CTA', 'CTT', 'TTA', 'TTG'],
	'AGC': ['AGT', 'TCG', 'TCA', 'TCC', 'TCT'],
	'ACA': ['ACC', 'ACG', 'ACT'],
	'AGA': ['AGG', 'CGA', 'CGC', 'CGG', 'CGT'],
	'CAT': ['CAC'],
	'AAT': ['AAC'],
	'ATT': ['ATC', 'ATA'],
	'CTG': ['CTA', 'CTC', 'CTT', 'TTA', 'TTG'],
	'CTA': ['CTG', 'CTC', 'CTT', 'TTA', 'TTG'],
	'ACT': ['ACC', 'ACA', 'ACG'],
	'CAC': ['CAT'],
	'ACG': ['ACC', 'ACA', 'ACT'],
	'CAA': ['CAG'],
	'AGT': ['AGC', 'TCG', 'TCA', 'TCC', 'TCT'],
	'CAG': ['CAA'],
	'CCG': ['CCT', 'CCA', 'CCC'],
	'CCC': ['CCT', 'CCG', 'CCA'],
	'TAT': ['TAC'],
	'GGT': ['GGG', 'GGA', 'GGC'],
	'TGT': ['TGC'],
	'CGA': ['AGG', 'AGA', 'CGC', 'CGG', 'CGT'],
	'CCA': ['CCT', 'CCG', 'CCC'],
	'TCT': ['AGC', 'AGT', 'TCG', 'TCA', 'TCC'],
	'GAT': ['GAC'],
	'CGG': ['AGG', 'AGA', 'CGA', 'CGC', 'CGT'],
	'TTT': ['TTC'],
	'TGC': ['TGT'],
	'GGG': ['GGT', 'GGA', 'GGC'],
	'TAG': ['TAA', 'TGA'],
	'GGA': ['GGT', 'GGG', 'GGC'],
	'TAA': ['TAG', 'TGA'],
	'GGC': ['GGT', 'GGG', 'GGA'],
	'TAC': ['TAT'],
	'GAG': ['GAA'],
	'TCG': ['AGC', 'AGT', 'TCA', 'TCC', 'TCT'],
	'TTA': ['CTG', 'CTA', 'CTC', 'CTT', 'TTG'],
	'GAC': ['GAT'],
	'TCC': ['AGC', 'AGT', 'TCG', 'TCA', 'TCT'],
	'GAA': ['GAG'],
	'TCA': ['AGC', 'AGT', 'TCG', 'TCC', 'TCT'],
	'GCA': ['GCC', 'GCG', 'GCT'],
	'GTA': ['GTC', 'GTG', 'GTT'],
	'GCC': ['GCA', 'GCG', 'GCT'],
	'GTC': ['GTA', 'GTG', 'GTT'],
	'GCG': ['GCA', 'GCC', 'GCT'],
	'GTG': ['GTA', 'GTC', 'GTT'],
	'TTC': ['TTT'],
	'GTT': ['GTA', 'GTC', 'GTG'],
	'GCT': ['GCA', 'GCC', 'GCG'],
	'ACC': ['ACA', 'ACG', 'ACT'],
	'TGA': ['TAG', 'TAA'],
	'TTG': ['CTG', 'CTA', 'CTC', 'CTT', 'TTA'],
	'CGT': ['AGG', 'AGA', 'CGA', 'CGC', 'CGG'],
	'TGG': [],
	'CGC': ['AGG', 'AGA', 'CGA', 'CGG', 'CGT']}

codon_freqs = { #It was a big headache to get these all to sum to one due to some rounding in the source data. Some values were arbitrarily bumped up or down by max 0.01 to get them to sum to one
	"A": {'yeast': [0.29, 0.22, 0.11, 0.38], 'neutral': [0.25, 0.25, 0.25, 0.25], 'codons': ['gca', 'gcc', 'gcg', 'gct'], 'human': [0.11, 0.4, 0.26, 0.23], 'coli': [0.23, 0.26, 0.33, 0.18], 'strep': [0.06, 0.17, 0.39, 0.38]},
	"C": {'yeast': [0.37, 0.63], 'neutral': [0.5, 0.5], 'codons': ['tgc', 'tgt'], 'human': [0.45, 0.55], 'coli': [0.54, 0.46], 'strep': [0.67, 0.33]},
	"E": {'yeast': [0.71, 0.29], 'neutral': [0.5, 0.5], 'codons': ['gaa', 'gag'], 'human': [0.42, 0.58], 'coli': [0.68, 0.32], 'strep': [0.74, 0.26]},
	"D": {'yeast': [0.35, 0.65], 'neutral': [0.5, 0.5], 'codons': ['gac', 'gat'], 'human': [0.46, 0.54], 'coli': [0.37, 0.63], 'strep': [0.68, 0.32]},
	"G": {'yeast': [0.22, 0.19, 0.12, 0.47], 'neutral': [0.25, 0.25, 0.25, 0.25], 'codons': ['gga', 'ggc', 'ggg', 'ggt'], 'human': [0.25, 0.34, 0.16, 0.25], 'coli': [0.13, 0.37, 0.15, 0.35], 'strep': [0.11, 0.17, 0.38, 0.34]},
	"F": {'yeast': [0.41, 0.59], 'neutral': [0.5, 0.5], 'codons': ['ttc', 'ttt'], 'human': [0.45, 0.55], 'coli': [0.42, 0.58], 'strep': [0.74, 0.26]},
	"I": {'yeast': [0.28, 0.26, 0.46], 'neutral': [0.3333333333333333, 0.3333333333333333, 0.3333333333333334], 'codons': ['ata', 'atc', 'att'], 'human': [0.36, 0.48, 0.16], 'coli': [0.12, 0.39, 0.49], 'strep': [0.6, 0.25, 0.15]},
	"H": {'yeast': [0.36, 0.64], 'neutral': [0.5, 0.5], 'codons': ['cac', 'cat'], 'human': [0.41, 0.59], 'coli': [0.43, 0.57], 'strep': [0.71, 0.29]},
	"K": {'yeast': [0.58, 0.42], 'neutral': [0.5, 0.5], 'codons': ['aaa', 'aag'], 'human': [0.42, 0.58], 'coli': [0.74, 0.26], 'strep': [0.78, 0.22]},
	"*": {'yeast': [0.47, 0.24, 0.29], 'neutral': [0.3333333333333333, 0.3333333333333333, 0.3333333333333334], 'codons': ['taa', 'tag', 'tga'], 'human': [0.28, 0.2, 0.52], 'coli': [0.61, 0.09, 0.3], 'strep': [0.74, 0.15, 0.11]},
	"M": {'yeast': [1.0], 'neutral': [1.0], 'codons': ['atg'], 'human': [1.0], 'coli': [1.0], 'strep': [1.0]},
	"L": {'yeast': [0.13, 0.06, 0.11, 0.13, 0.28, 0.29], 'neutral': [0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666667], 'codons': ['cta', 'ctc', 'ctg', 'ctt', 'tta', 'ttg'], 'human': [0.06, 0.13, 0.13, 0.41, 0.07, 0.2], 'coli': [0.04, 0.1, 0.47, 0.12, 0.14, 0.13], 'strep': [0.33, 0.19, 0.19, 0.06, 0.17, 0.06]},
	"N": {'yeast': [0.41, 0.59], 'neutral': [0.5, 0.5], 'codons': ['aac', 'aat'], 'human': [0.46, 0.54], 'coli': [0.51, 0.49], 'strep': [0.64, 0.36]},
	"Q": {'yeast': [0.69, 0.31], 'neutral': [0.5, 0.5], 'codons': ['caa', 'cag'], 'human': [0.25, 0.75], 'coli': [0.34, 0.66], 'strep': [0.74, 0.26]},
	"P": {'yeast': [0.42, 0.15, 0.12, 0.31], 'neutral': [0.25, 0.25, 0.25, 0.25], 'codons': ['cca', 'ccc', 'ccg', 'cct'], 'human': [0.12, 0.33, 0.28, 0.27], 'coli': [0.2, 0.13, 0.49, 0.18], 'strep': [0.06, 0.06, 0.43, 0.45]},
	"S": {'yeast': [0.11, 0.16, 0.21, 0.16, 0.1, 0.26], 'neutral': [0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666667], 'codons': ['agc', 'agt', 'tca', 'tcc', 'tcg', 'tct'], 'human': [0.18, 0.06, 0.15, 0.15, 0.22, 0.24], 'coli': [0.24, 0.16, 0.14, 0.15, 0.14, 0.17], 'strep': [0.31, 0.03, 0.3, 0.19, 0.06, 0.11]},
	"R": {'yeast': [0.47, 0.21, 0.07, 0.06, 0.04, 0.15], 'neutral': [0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666666, 0.16666666666666667], 'codons': ['aga', 'agg', 'cga', 'cgc', 'cgg', 'cgt'], 'human': [0.22, 0.08, 0.11, 0.2, 0.19, 0.2], 'coli': [0.06, 0.04, 0.07, 0.36, 0.11, 0.36], 'strep': [0.06, 0.37, 0.1, 0.05, 0.12, 0.3]},
	"T": {'yeast': [0.3, 0.22, 0.13, 0.35], 'neutral': [0.25, 0.25, 0.25, 0.25], 'codons': ['aca', 'acc', 'acg', 'act'], 'human': [0.12, 0.36, 0.24, 0.28], 'coli': [0.16, 0.4, 0.25, 0.19], 'strep': [0.15, 0.14, 0.3, 0.41]},
	"W": {'yeast': [1.0], 'neutral': [1.0], 'codons': ['tgg'], 'human': [1.0], 'coli': [1.0], 'strep': [1.0]},
	"V": {'yeast': [0.21, 0.21, 0.19, 0.39], 'neutral': [0.25, 0.25, 0.25, 0.25], 'codons': ['gta', 'gtc', 'gtg', 'gtt'], 'human': [0.47, 0.24, 0.18, 0.11], 'coli': [0.17, 0.2, 0.35, 0.28], 'strep': [0.21, 0.16, 0.4, 0.23]},
	"Y": {'yeast': [0.44, 0.56], 'neutral': [0.5, 0.5], 'codons': ['tac', 'tat'], 'human': [0.43, 0.57], 'coli': [0.41, 0.59], 'strep': [0.71, 0.29]}}

neutral_transformation_matrix = {
	'ACC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
	'ATG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'AAG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'AAA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'ATC': {'ACC': 1.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'AAC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 1.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'ATA': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'AGG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
	'CCT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 1.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
	'CTC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'AGC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'ACA': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'AGA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CAT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'AAT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'ATT': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CTG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CTA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'ACT': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
	'CAC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'ACG': {'ACC': 1.0, 'ATG': 1.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CAA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'AGT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 1.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CCA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CCG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CCC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
	'TAT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
	'GGT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TGT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 1.0},
	'CGA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CAG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CGC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'GAT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 1.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CGG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 1.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
	'CTT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 1.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 1.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TGC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 1.0, 'TCT': 0.0},
	'GGG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
	'TAG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 1.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
	'GGA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TAA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'GGC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 1.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TAC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
	'GAG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 1.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 1.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TCG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 1.0, 'TGG': 1.0, 'TCT': 1.0},
	'TTA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TTT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
	'GAC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 1.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 1.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 1.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'CGT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 1.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 1.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 1.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 1.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'GAA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 1.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 1.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TCA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 1.0},
	'GCA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 1.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 1.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'GTA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 1.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 1.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 1.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'GCC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
	'GTC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 1.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 1.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 1.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'GCG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 1.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 1.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'GTG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 1.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 1.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TTC': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 1.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 1.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 1.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0},
	'GTT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 1.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 1.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 1.0, 'GCC': 0.0, 'GTC': 1.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'GCT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 1.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 1.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 1.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 1.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 1.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
	'TGA': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 1.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 1.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 1.0, 'TAA': 1.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 0.0, 'TTA': 1.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
	'TTG': {'ACC': 0.0, 'ATG': 1.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 1.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 1.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 1.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 1.0, 'TCT': 0.0},
	'TCC': {'ACC': 1.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 1.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 0.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 1.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 1.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 1.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 1.0},
	'TGG': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 1.0, 'CCT': 0.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 0.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 0.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 1.0, 'CTT': 0.0, 'TGC': 1.0, 'GGG': 1.0, 'TAG': 1.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 0.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 0.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 0.0, 'TGA': 1.0, 'TTG': 1.0, 'TCC': 0.0, 'TGG': 0.0, 'TCT': 0.0},
	'TCT': {'ACC': 0.0, 'ATG': 0.0, 'AAG': 0.0, 'AAA': 0.0, 'ATC': 0.0, 'AAC': 0.0, 'ATA': 0.0, 'AGG': 0.0, 'CCT': 1.0, 'CTC': 0.0, 'AGC': 0.0, 'ACA': 0.0, 'AGA': 0.0, 'CAT': 0.0, 'AAT': 0.0, 'ATT': 0.0, 'CTG': 0.0, 'CTA': 0.0, 'ACT': 1.0, 'CAC': 0.0, 'ACG': 0.0, 'CAA': 0.0, 'AGT': 0.0, 'CCA': 0.0, 'CCG': 0.0, 'CCC': 0.0, 'TAT': 1.0, 'GGT': 0.0, 'TGT': 1.0, 'CGA': 0.0, 'CAG': 0.0, 'CGC': 0.0, 'GAT': 0.0, 'CGG': 0.0, 'CTT': 0.0, 'TGC': 0.0, 'GGG': 0.0, 'TAG': 0.0, 'GGA': 0.0, 'TAA': 0.0, 'GGC': 0.0, 'TAC': 0.0, 'GAG': 0.0, 'TCG': 1.0, 'TTA': 0.0, 'TTT': 1.0, 'GAC': 0.0, 'CGT': 0.0, 'GAA': 0.0, 'TCA': 1.0, 'GCA': 0.0, 'GTA': 0.0, 'GCC': 0.0, 'GTC': 0.0, 'GCG': 0.0, 'GTG': 0.0, 'TTC': 0.0, 'GTT': 0.0, 'GCT': 1.0, 'TGA': 0.0, 'TTG': 0.0, 'TCC': 1.0, 'TGG': 0.0, 'TCT': 0.0}}

if __name__ == "__main__":
	
	from optparse import OptionParser
	import sys, bisect, random, os, numpy

	parser = OptionParser()
	parser.add_option('--setup', action = 'store', type = 'string', dest = 'setup_file', default = './setup.csv', help = "create a setup file from the template to put here")
	(option, args) = parser.parse_args()
	
	main(option.setup_file)