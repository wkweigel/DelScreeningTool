#These should be the specific primers used for a specif
# Dictionary of first PCR primer. These should be found at the beginning
# (5' end) of each read.
pcr_primer_1 = {
	'ACAGCA': '1a-02',
	'ACATGT': '1a-03',
	'ACGACG': '1a-04',
	'ACGCGA': '1a-05',
	'ACTAGC': '1a-06',
	'ACTCTG': '1a-07',
	'ACTGAT': '1a-08',
	'AGACTA': '1a-09',
	'AGAGAG': '1a-10',
    'TAGATG': '1a-11',
	'TAGTGT': '1a-12',
	'TATCAC': '1a-13',
}

# Dictionary of second PCR primer. These should be found at the end
# (3' end) of each read. Detection may be difficult due to indels.
pcr_primer_2 = {
	'ATATCG': '1b-a',
	'ATCAGT': '1b-b',
	'ATGATA': '1b-c',
	'TGTCAG': '1b-d',
	'TTGTGC': '1b-e',
}
