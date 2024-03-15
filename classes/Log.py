import datetime

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
	
	def generate_log(self, log_path):
		with open(log_path, "w+") as log_file:
			heading1='User Inputs'
			log_file.write(f'{heading1:=^50}\n')
			for k,v in self.user_inputs.items():
					log_file.write(f'{k:<20} {str(v)}\n')
			log_file.write(f'\n')

			heading2='Data Pre-Processing'
			log_file.write(f'{heading2:=^50}\n')
			for k,v in self.preprocessing_dict.items():
				log_file.write(f'{k:<20} {str(v)}\n')
			log_file.write(f'\n')

			heading3='Code Correction'
			log_file.write(f'{heading3:=^50}\n')

			log_file.write(f'Rows before correction: {self.rows_before} \n')
			for code, info in self.correction_info.items():
				log_file.write(f'{code:<4}| Corrected: {info[0]:<7} Rejected: {info[1]:<7} Attrition Rate: {info[2]}\n')
			log_file.write(f'Rows after correction: {self.rows_after} \n')
			log_file.write('\n')
			log_file.write('Correction times: \n')
			for code, time in self.correction_times.items():
				log_file.write(f'{code:<5} | {time} sec \n')


