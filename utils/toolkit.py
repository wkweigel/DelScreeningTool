from Bio import pairwise2
import datetime
import numpy as np
import pandas as pd


def perform_alignment_with_score(main_sequence=str, subsequence=str):
    """Perform an alignment for two string sequences and return the score"""
    #Perform alignment using +1 for matches and -1 for mismatches, insertions, or deletions
    alignments = pairwise2.align.localms(main_sequence, subsequence, 1, -1, -1, -1, one_alignment_only=True)
    if alignments:
        best_alignment = alignments[0]
        score=best_alignment[-3]
        return(score)


def perform_alignment_with_report(long_sequence=str, subsequence=str, QUERY_TOLERANCE=int):
	"""Perform an alignment for two string sequences and return five reporter values:
		start position, end position, insertions, deletions, and mismatches"""
	#Perform alignment using +1 for matches and -1 for mismatches, insertions, or deletions
	alignments = pairwise2.align.localms(long_sequence, subsequence, 1, -1, -1, -1, one_alignment_only=True)

	#Define the alignments for each sequence
	#Note: "-" in the long_sequence indicates an deletion, "-" in the subsequence indicates an insertion
	if alignments:
		best_alignment = alignments[0]
		aligned_long_sequence = best_alignment[0]
		aligned_subsequence = best_alignment[1]

		#Extract the start and endpoints of the alignment on the long_sequence
		start = best_alignment[-2]
		end = best_alignment[-1]

		#Determine the number and type of errors in the alignment
		insertions=0
		deletions=0
		mismatches=0

		for i in range(len(aligned_subsequence[start:end])):
			b1=aligned_long_sequence[start:end][i]
			b2=aligned_subsequence[start:end][i]
			if b1 == "-":
				deletions+=1
			if b2 == "-":
				insertions+=1
			if b1 != b2 and b1 != '-' and b2 != '-':
				mismatches+=1

		return(start, end-deletions, insertions, deletions, mismatches)


def extract_code_from_sequence(sequence=str, query_seq=str, code_length=int, QUERY_TOLERANCE=int):
	"""Locate and extract an encoding sequence from a parent sequence based on a query.
	
	Arguments
	==========
	sequence: A sequence containing an encoding subsequence.
	query_seq: A static sequence immediately preceding the desired encoding sequence.
	code_length: The length of the encoding sequence.

	Returns
	=======
	A single encoding sequence.
	
	"""

	s,e,i,d,m=perform_alignment_with_report(sequence, query_seq, QUERY_TOLERANCE)
	if i+d+m<=QUERY_TOLERANCE: #tolerance of 3 errors on alignment
		code=sequence[e:e+code_length]
		code=(code)
	else:
		code="No Match"
	return(code)


# From a list of sequences, return a dict of the unique sequences counts
def get_seq_counts(SeqList=list):
	"""Compile a dict of the counts for each unique sequence in a list."""
	count_dict = {}
	for i in SeqList:
		count_dict[i] = count_dict.get(i, 0) + 1
	return(count_dict)


def count_matching_characters(str1, str2):
	"""Counts the number of matching characters between two strings of equal length."""
	# Use a list comprehension to create a list of 0s and 1s 
	matching_characters = [1 if char1 == char2 else 0 for char1, char2 in zip(str1, str2)]
	# Use the sum function to count the matching characters
	count = sum(matching_characters)
	return count


def get_closest_match(incorrect_sequence, correct_sequence_list):
	"""For an incorrect sequence, find the closest matching sequence from a list of correct sequences.
	
	Arguments
	==========
	incorrect_sequence: An encoding sequence with at least one error
	correct_sequence_list: A list of the possible correct encoding sequences
	
	Returns
	=======
	num_errors: An int value representing the total number of errors in the sequence.
	best_match: A string corresponding to the closest matching correct sequence.
	"""
	score_dict={}
	if incorrect_sequence in correct_sequence_list: #If the sequence is correct
		num_matches=len(incorrect_sequence)
		best_match=incorrect_sequence
	else:
		for correct_sequence in correct_sequence_list:
			score=count_matching_characters(correct_sequence, incorrect_sequence)
			score_dict[correct_sequence]=score
		best_score=max(score_dict.values())
		best_match=max(score_dict, key=lambda k: score_dict.get(k))
		num_matches=best_score
	num_errors=len(incorrect_sequence)-num_matches
	
	return(num_errors, best_match)


def calc_variance(df):
	"""Calculate the multiplicative variance using the counts of the unique values in column of a dataframe.
	  
	Arguments
	==========
	df: The dataframe used to calculate the column-wise variance for.

	Returns
	=======
	variance_scores: A list of scores for each column in the input dataframe.
		
	"""
	variance_scores=[]
	for i in range(len(df.columns)):
		variance=1
		value_counts = df[i+1].value_counts()
		for val in value_counts:
			variance=variance*val
		variance_scores.append(variance)

	return(variance_scores)


def update_code_df(code_df:pd.DataFrame, code_name:str, code_list:list):
	code_df[[f'{code_name}_errors', f'corrected_{code_name}']] = code_df[f'{code_name}'].progress_apply(lambda code: pd.Series(get_closest_match(code, code_list)))
	return(code_df)



def filter_code_df(code_df:pd.DataFrame, corrected_df:pd.DataFrame ,code_name:str, correction_dict:dict):
	corrected_df_before_rows=corrected_df.shape[0] #rows before starting current round of filtering 
	corrected_df_initial_rows=code_df.shape[0] #rows before starting any filtering 
	
	corrected=np.count_nonzero(code_df[f'{code_name}_errors']==1) #count the rows that were corrected
	rejected=np.count_nonzero(code_df[f'{code_name}_errors']>1) #count the rows that were rejected

	current_corrected=np.count_nonzero(corrected_df[f'{code_name}_errors']==1) #count the rows that currently corrected
	current_rejected=np.count_nonzero(corrected_df[f'{code_name}_errors']>1) #count the rows that currently rejected

	corrected_df=corrected_df[corrected_df[f'{code_name}_errors'] <=1] #filter out any rows with too many errors
	
	corrected_df_remaining_rows=corrected_df.shape[0]

	#attrition=round(corrected_df_remaining_rows/corrected_df_initial_rows, ndigits=2)
	attrition=round(((corrected_df_initial_rows-((corrected_df_initial_rows-corrected_df_before_rows)+current_rejected))/corrected_df_initial_rows), ndigits=2)
	correction_dict[code_name]=(corrected, rejected, attrition)
	return((corrected_df, correction_dict))
	

class Logging:
	def __init__(self):
		#self.output_path=output_path
		self.init_time=datetime.datetime.now()
		self.user_inputs={}
		
	def store_preprocessing_records(self, n_before:int, n_after:int, time:float):
			self.n_reads_before=n_before
			self.n_reads_after=n_after
			self.n_reads_removed=n_before-n_after
			self.preprocessing_dict={'Total reads': n_before, 
									 'Extracted reads':f'{n_after} ({round((n_after/n_before)*100, ndigits=1)}%)',
									 'Removed reads' : f'{self.n_reads_removed} ({round((self.n_reads_removed/n_before)*100, ndigits=1)}%)',
									 'Preproc time' : f'{time} sec'
									 }
		
	def store_extraction_records(self, start_time:float, end_time:float):
			
			self.extraction_duration=round(start_time-end_time, ndigits=4)
		
	def store_correction_records(self, n_before:int, n_after:int, info_dict:dict, time_dict:dict):
			self.rows_before=n_before
			self.rows_after=n_after
			self.correction_info=info_dict
			self.correction_times=time_dict

			print(f'Rows before correction: ', n_before)
			for code, info in info_dict.items():
				print(f'{code}| Corrected: {info[0]}, Rejected: {info[1]}, Attrition Rate: {info[2]}')
			print(f'Rows after correction: ', n_after)



class Codes:
	def __init__(self, bb_codes:dict, pcr_codes:dict):
		self.bb_codes = bb_codes
		self.pcr_codes = pcr_codes

		#Create lists of the correct encoding sequences used for error checking
		#These are defined with the keys from dicts found in the nadel_fragments.py file and the primers.py file
		code_sequences_bb1=list(fragment_set_1.keys())
		code_sequences_bb2=list(fragment_set_2.keys())
		code_sequences_pcr1=list(pcr_primer_1.keys())
		code_sequences_pcr2=list(pcr_primer_2.keys())