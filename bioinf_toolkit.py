from bio_structs import *


def validate_seq(seq):
	seq = seq.upper()
	for nucleotide in seq:
		if nucleotide not in Nucleotides['DNA']:
			return False
		else:
			return True


def count_nuc_frequency(seq):
	seq = seq.upper()
	tempFreq = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0}
	for nucleotide in seq:
		tempFreq[nucleotide] += 1
	return tempFreq


def transcription(seq):
	return seq.replace('T', 'U')


def complement(seq):
	''' Takes a DNA strand string and finds its opposite base pair match. '''
	table = ''.maketrans('ACTG', 'TGAC')
	return seq.translate(table)


def translation(seq, init_index=0):
	# i - position in the gene, init_index - start nucleotide
	return [DNA_Codons[seq[i:i+3]] for i in range(init_index, len(seq) - 2, 3)]


def gc_content(seq):
	return (seq.count('G') + seq.count('C')) / len(seq)


def gc_content_of_subseq(seq, k=20):
	''' Find GC content in a DNA/RNA of all subsequences lenght k. '''
	result = []
	for i in range(0, len(seq) - k + 1, k):
		subseq = seq[i:i + k]
		result.append(gc_content(subseq))
	return result

def codon_freq(seq):
	