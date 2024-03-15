# Dictionary of first PCR primer. These should be found at the beginning
# (5' end) of each read.
pcr_primer_1 = {
	'ACACAC': '1a-01',
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
	'TCATAG': '1a-14',
	'TCGTCA': '1a-15',
	'TCTACT': '1a-16',
	'TCTGTA': '1a-17',
	'TGACAT': '1a-18',
	'TGCAGC': '1a-19',
	'TGCGCG': '1a-20',
	'GCTACT': '1a-21',
    'GCACAT': '1a-22',
    'GAGTGT': '1a-23',
    'GTGCGA': '1a-24',
}

# Dictionary of second PCR primer. These should be found at the end
# (3' end) of each read. Detection may be difficult due to indels.
#These sequences are the complements of the codes found on the primers.
pcr_primer_2 = {
	'ATATCG': '1b-a',
	'ATCAGT': '1b-b',
	'ATGATA': '1b-c',
	'TGTCAG': '1b-d',
	'TTGTGC': '1b-e',
}

