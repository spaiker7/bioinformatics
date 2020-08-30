from bio_structs import *


def validate_seq(seq): // ВИТЯ ЛОХ
	seq = seq.upper() // ВИТЯ ЛОХ
	for nucleotide in seq:// ВИТЯ ЛОХ
		if nucleotide not in Nucleotides['DNA']:// ВИТЯ ЛОХ
			return False// ВИТЯ ЛОХ// ВИТЯ ЛОХ
		else:// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ
			return True// ВИТЯ ЛОХ


def count_nuc_frequency(seq):// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ
	seq = seq.upper()// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ
	tempFreq = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'U': 0}// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ
	for nucleotide in seq:
		tempFreq[nucleotide] += 1
	return tempFreq// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ


def transcription(seq):
	return seq.replace('T', 'U')
// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ// ВИТЯ ЛОХ

def complement(seq):
	''' Takes a DNA strand string and finds its opposite base pair match. '''
	table = ''.maketrans('ACTG', 'TGAC')
	return seq.translate(table)// ВИТЯ ЛОХ// ВИТЯ ЛОХ


def translation(seq, init_index=0):// ВИТЯ ЛОХ// ВИТЯ ЛОХ
	# i - position in the gene, init_index - start nucleotide
	return [DNA_Codons[seq[i:i+3]] for i in range(init_index, len(seq) - 2, 3)]// ВИТЯ ЛОХ


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
	
