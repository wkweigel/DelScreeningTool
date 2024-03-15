from Bio import SeqIO
from tqdm import tqdm
import pandas as pd
import os
from .context import library

#Import encoding sequence dictionaries from the fragment_codes.py file
from library.fragment_codes import code_set_1, code_set_2, code_set_3

#Import BB smiles dictionaries from the fragment_smiles.py file
from library.fragment_smiles import smiles_set_1, smiles_set_2, smiles_set_3

#Import  encoding sequence dictionaries from the primers_screen.py file
from library.primers_screen import pcr_primer_1, pcr_primer_2


class Codes:
    '''A class for defining, extracting, and tabulating encoding sequences'''
    def __init__(self, n_cycles:int):
        # The number of library cycles (building blocks) to consider 
        self.n_cycles=n_cycles

        # A dict used to hold the sequences for each encoding region
        self.code_dict={}

        # Library uses 2 primers regardless of cycle number
        self.pcr1_code_dict = pcr_primer_1
        self.code_dict['pcr1']=list(pcr_primer_1.keys())
        self.pcr2_code_dict = pcr_primer_2
        self.code_dict['pcr2']=list(pcr_primer_2.keys())

        # Create class objects based on the number of cycles (n_cycles) in the library
        if n_cycles>=1:
            self.bb1_code_dict = code_set_1
            self.bb1_smiles_dict = smiles_set_1
            self.code_dict['bb1']=list(code_set_1.keys())
            
        if n_cycles>=2:
            self.bb2_code_dict = code_set_2
            self.bb2_smiles_dict = smiles_set_2
            self.code_dict['bb2']=list(code_set_2.keys())
            
        if n_cycles>=3:
            self.bb3_code_dict = code_set_3
            self.bb3_smiles_dict = smiles_set_3
            self.code_dict['bb3']=list(code_set_3.keys())
        
        self.user_inputs=None
    
    def load_user_inputs(self, user_inputs:dict):
        self.user_inputs=user_inputs
        subfolder=self.user_inputs['NAME']
        self.output_location=f'outputs/{subfolder}'
        os.makedirs(self.output_location, exist_ok=True)
    
    def preprocess_fastq(self):
        if self.user_inputs==None:
            print(f'User inputs have not been loaded. Please run .load_user_inputs')    
        else:
            fastq_file=self.user_inputs['RAW_FASTQ_FILE']
            filename=self.user_inputs['PROC_FASTQ_FILE']
            with open(f'{self.output_location}/{filename}', "w+") as output_file:
                n_records=0
                w_records=0
                # Use Biopython to iterate over the records in the FASTQ file
                for record in tqdm(SeqIO.parse(fastq_file, "fastq")):
                    n_records+=1
                    if self.user_inputs['ADAPT_SEQ'] in str(record.seq):
                        output_file.write(str(record.seq) + "\n")
                        w_records+=1
                        if self.user_inputs['TEST_RUN']==True and w_records>=self.user_inputs['N_TEST_SEQUENCES']:
                            break
            return(n_records, w_records)
        
    def extract_codes(self, generate_csv=True, csv_filename='Extracted_Codes.csv'):
        '''Extract library encoding sequences from a preprocessed .txt file.'''
        pcr1_codes=[]
        pcr2_codes=[]
        bb1_codes=[]
        bb2_codes=[]
        bb3_codes=[]

        bb_code_len=self.user_inputs['BB_ENCODING_LEN']
        pcr_code_len=self.user_inputs['PCR_ENCODING_LEN']
        filename=self.user_inputs['PROC_FASTQ_FILE']
        with open(f'{self.output_location}/{filename}', 'r') as seq_file:
                line_count = len(seq_file.readlines())
                print(f"Processing {line_count} sequences...")

        with open(f'{self.output_location}/{filename}', 'r') as seq_file:
            for record in tqdm(seq_file):
                pcr1_codes.append(record[:self.user_inputs['PCR_ENCODING_LEN']]) #No start idx is needed if the pcr1 is the first n bases in the sequence read 
                pcr2_codes.append(record[self.user_inputs['PCR2_START_IDX'] : self.user_inputs['PCR2_START_IDX']  + pcr_code_len])
                if self.n_cycles>=1:
                    bb1_codes.append(record[self.user_inputs['BB1_START_IDX'] : self.user_inputs['BB1_START_IDX'] + bb_code_len])
                if self.n_cycles>=2:
                    bb2_codes.append(record[self.user_inputs['BB2_START_IDX'] : self.user_inputs['BB2_START_IDX'] + bb_code_len])
                if self.n_cycles>=3:
                    bb3_codes.append(record[self.user_inputs['BB3_START_IDX'] : self.user_inputs['BB3_START_IDX'] + bb_code_len])
        
        if self.n_cycles==1:
            self.code_df = pd.DataFrame(list(zip(bb1_codes, pcr1_codes, pcr2_codes)),columns=['bb1', 'pcr1', 'pcr2'])
        if self.n_cycles==2:
            self.code_df = pd.DataFrame(list(zip(bb1_codes, bb2_codes, pcr1_codes, pcr2_codes)),columns=['bb1', 'bb2', 'pcr1', 'pcr2'])
        if self.n_cycles==3:
            self.code_df = pd.DataFrame(list(zip(bb1_codes, bb2_codes, bb3_codes, pcr1_codes, pcr2_codes)),columns=['bb1', 'bb2', 'bb3', 'pcr1', 'pcr2'])
        
        if generate_csv == True:
            self.code_df.to_csv(f"{self.output_location}/{csv_filename}")
        
        



        