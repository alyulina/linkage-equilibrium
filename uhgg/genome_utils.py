import os
import math
import pandas as pd
import numpy as np

complement_nuc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
def reverse_complement(seq):
    return ''.join(complement_nuc[i] for i in seq)[::-1]

codon_table = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R', 'AAT':'N', 'AAC':'N', 'GAT':'D', 'GAC':'D', 'TGT':'C', 'TGC':'D', 'CAA':'Q', 'CAG':'Q', 'GAA':'E', 'GAG':'E', 'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 'CAT':'H', 'CAC':'H', 'ATT':'I', 'ATC':'I', 'ATA':'I', 'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 'AAA':'K', 'AAG':'K', 'ATG':'M', 'TTT':'F', 'TTC':'F', 'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S', 'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 'TGG':'W', 'TAT':'Y', 'TAC':'Y', 'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V', 'TAA':'!', 'TGA':'!', 'TAG':'!'}
codon_degeneracy_table = {'GCT': [1, 1, 4], 'GCC': [1, 1, 4], 'GCA': [1, 1, 4], 'GCG': [1, 1, 4], 'CGT': [1, 1, 4], 'CGC': [1, 1, 4], 'CGA': [2, 1, 4], 'CGG': [2, 1, 4], 'AGA': [2, 1, 2], 'AGG': [2, 1, 2], 'AAT': [1, 1, 2], 'AAC': [1, 1, 2], 'GAT': [1, 1, 2], 'GAC': [1, 1, 2], 'TGT': [1, 1, 1], 'TGC': [1, 1, 1], 'CAA': [1, 1, 2], 'CAG': [1, 1, 2], 'GAA': [1, 1, 2], 'GAG': [1, 1, 2], 'GGT': [1, 1, 4], 'GGC': [1, 1, 4], 'GGA': [1, 1, 4], 'GGG': [1, 1, 4], 'CAT': [1, 1, 2], 'CAC': [1, 1, 2], 'ATT': [1, 1, 3], 'ATC': [1, 1, 3], 'ATA': [1, 1, 3], 'TTA': [2, 1, 2], 'TTG': [2, 1, 2], 'CTT': [1, 1, 4], 'CTC': [1, 1, 4], 'CTA': [2, 1, 4], 'CTG': [2, 1, 4], 'AAA': [1, 1, 2], 'AAG': [1, 1, 2], 'ATG': [1, 1, 1], 'TTT': [1, 1, 2], 'TTC': [1, 1, 2], 'CCT': [1, 1, 4], 'CCC': [1, 1, 4], 'CCA': [1, 1, 4], 'CCG': [1, 1, 4], 'TCT': [1, 1, 4], 'TCC': [1, 1, 4], 'TCA': [1, 1, 4], 'TCG': [1, 1, 4], 'AGT': [1, 1, 2], 'AGC': [1, 1, 2], 'ACT': [1, 1, 4], 'ACC': [1, 1, 4], 'ACA': [1, 1, 4], 'ACG': [1, 1, 4], 'TGG': [1, 1, 1], 'TAT': [1, 1, 2], 'TAC': [1, 1, 2], 'GTT': [1, 1, 4], 'GTC': [1, 1, 4], 'GTA': [1, 1, 4], 'GTG': [1, 1, 4], 'TAA': [1, 2, 2], 'TGA': [1, 2, 1], 'TAG': [1, 1, 2]}
# THIS IS BEN'S CODE
base_table = {'A':'T','T':'A','G':'C','C':'G'}
codon_degeneracy_table = {}
codon_annotation_table = {}
codon_substitution_table = {}

# calculate opportunities for each codon
for codon in codon_table:
    codon_annotation_table[codon] = {}
    codon_substitution_table[codon] = {}

    fold_degeneracy = [0, 0, 0]
    for i in range(0, 3):
        codon_annotation_table[codon][i] = {}
        codon_substitution_table[codon][i] = {}

        codon_as_list = list(codon)
        for base in base_table:
            codon_as_list[i] = base
            new_codon = "".join(codon_as_list)
            codon_substitution_table[codon][i][base] = new_codon
            # print codon, new_codon
            if codon_table[codon] == codon_table[new_codon]:
                # synonymous change!
                fold_degeneracy[i] += 1
                codon_annotation_table[codon][i][base] = 'S'
            elif codon_table[codon] == '!' or codon_table[new_codon] == '!':  # preumature stop of loss of stop
                codon_annotation_table[codon][i][base] = 'N'
            else:
                codon_annotation_table[codon][i][base] = 'M'

    codon_degeneracy_table[codon] = fold_degeneracy

print(codon_degeneracy_table)
print(codon_annotation_table)
print(codon_substitution_table)



def get_genome_sequence():
    return

def get_genome_annotation():
    return

def get_gene():
    return

def get_gene_sequence():
    return

