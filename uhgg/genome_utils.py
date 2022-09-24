import numpy as np
import pandas as pd

complement_nuc = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
def reverse_complement(seq):
    """Returns the reverse complement of a nucleotide sequence.

    Args:
        seq: nucleotide sequence in uppercase.
             type: str

    Returns: reverse complement nucleotide sequence in uppercase.
             type: str
    """
    # to add later: make it work for the lowercase as well, raise an exception if nucleotide is not in complement_nuc
    return ''.join(complement_nuc[x] for x in seq)[::-1]


codon_table = {'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
               'AGA': 'R', 'AGG': 'R', 'AAT': 'N', 'AAC': 'N', 'GAT': 'D', 'GAC': 'D', 'TGT': 'C', 'TGC': 'D',
               'CAA': 'Q', 'CAG': 'Q', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
               'CAT': 'H', 'CAC': 'H', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L',
               'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'AAA': 'K', 'AAG': 'K', 'ATG': 'M', 'TTT': 'F', 'TTC': 'F',
               'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
               'AGT': 'S', 'AGC': 'S', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'TGG': 'W', 'TAT': 'Y',
               'TAC': 'Y', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'TAA': '!', 'TGA': '!', 'TAG': '!'}

codon_degeneracy_table = {'GCT': [1, 1, 4], 'GCC': [1, 1, 4], 'GCA': [1, 1, 4], 'GCG': [1, 1, 4],
                          'CGT': [1, 1, 4], 'CGC': [1, 1, 4], 'CGA': [2, 1, 4], 'CGG': [2, 1, 4],
                          'AGA': [2, 1, 2], 'AGG': [2, 1, 2], 'AAT': [1, 1, 2], 'AAC': [1, 1, 2],
                          'GAT': [1, 1, 2], 'GAC': [1, 1, 2], 'TGT': [1, 1, 1], 'TGC': [1, 1, 1],
                          'CAA': [1, 1, 2], 'CAG': [1, 1, 2], 'GAA': [1, 1, 2], 'GAG': [1, 1, 2],
                          'GGT': [1, 1, 4], 'GGC': [1, 1, 4], 'GGA': [1, 1, 4], 'GGG': [1, 1, 4],
                          'CAT': [1, 1, 2], 'CAC': [1, 1, 2], 'ATT': [1, 1, 3], 'ATC': [1, 1, 3],
                          'ATA': [1, 1, 3], 'TTA': [2, 1, 2], 'TTG': [2, 1, 2], 'CTT': [1, 1, 4],
                          'CTC': [1, 1, 4], 'CTA': [2, 1, 4], 'CTG': [2, 1, 4], 'AAA': [1, 1, 2],
                          'AAG': [1, 1, 2], 'ATG': [1, 1, 1], 'TTT': [1, 1, 2], 'TTC': [1, 1, 2],
                          'CCT': [1, 1, 4], 'CCC': [1, 1, 4], 'CCA': [1, 1, 4], 'CCG': [1, 1, 4],
                          'TCT': [1, 1, 4], 'TCC': [1, 1, 4], 'TCA': [1, 1, 4], 'TCG': [1, 1, 4],
                          'AGT': [1, 1, 2], 'AGC': [1, 1, 2], 'ACT': [1, 1, 4], 'ACC': [1, 1, 4],
                          'ACA': [1, 1, 4], 'ACG': [1, 1, 4], 'TGG': [1, 1, 1], 'TAT': [1, 1, 2],
                          'TAC': [1, 1, 2], 'GTT': [1, 1, 4], 'GTC': [1, 1, 4], 'GTA': [1, 1, 4],
                          'GTG': [1, 1, 4], 'TAA': [1, 2, 2], 'TGA': [1, 2, 1], 'TAG': [1, 1, 2]}


def get_genome_sequence(fna_path):
    """Reads an .fna file and returns the genome sequence.

    Args:
        fna_path: path to the .fna file.
                  type: str

    Returns: a dictionary of nucleotide sequences as strings in uppercase with contig names as keys.
             type: {str: str}
    """

    with open(fna_path, 'r') as f:
        lines = f.readlines()

    contig_sequences = {}
    for i in lines:
        if i[0] == '>': # new contig
            contig_name = i.strip()[1:]
            if contig_name in contig_sequences.keys():
                raise ValueError('Duplicate contig!')
            else:
                contig_sequences[contig_name] = ''
        else:
            contig_sequences[contig_name] += i.strip()
    return contig_sequences


def get_genome_annotation(gff_path):
    """Reads a (preprocessed) .gff file and returns a pd.dataframe with genome annotation.

    Args:
        gff_path: path to the .gff file.
                  type: str

    Returns: a table containing the genome annotation.
             type: pd.dataframe

    Dependencies: pandas as pd
    """

    with open(gff_path, 'r') as f:
        lines = f.readlines()
    contigs = []
    sources = []
    types = []
    starts = []
    ends = []
    scores = []
    strands = []
    frames = []
    attributes = []

    for line in lines:
            items = line.split('\t')
            contigs.append(items[0])
            sources.append(items[1])
            types.append(items[2])
            starts.append(int(items[3]))
            ends.append(int(items[4]))
            scores.append(items[5])
            strands.append(items[6])
            frames.append(items[7])
            attributes.append(items[8])

    annotation = pd.DataFrame(list(zip(contigs, sources, types, starts, ends, scores, strands, frames, attributes)),
                      columns=['Contig', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute'])
    annotation['Gene ID'] = np.arange(annotation.shape[0])  # gene id assigned by order in the .gff file
    annotation['UHGG ID'] = annotation['Attribute'].apply(lambda x: x.split(';')[0].split('=')[-1])

    return annotation


def get_gene_location(gff_path, gene_id):
    """Gets gene coordinates, strand, and contig.

    Args:
        gff_path: path to the .gff genome annotation file.
                  type: str
        gene_id: gene id in the .gff file
                 type: int or str

    Returns: contig, gene start, end, and strand.
             type: str, int, int, str

    Dependencies: pandas as pd, get_genome_annotation(gff_path)
    """
    genome_annotation = get_genome_annotation(gff_path)
    gene_annotation = genome_annotation.loc[genome_annotation['Gene ID'] == gene_id]
    return gene_annotation.iloc[0, 0], gene_annotation.iloc[0, 3], gene_annotation.iloc[0, 4], gene_annotation.iloc[0, 6]


def get_gene_sequence(fna_path, gff_path, gene_id):
    """Gets the sequence of a gene from the reference genome.

    Args:
        fna_path: path to the .fna reference genome file.
                  type: str
        gff_path: path to the .gff file.
                  type: str
        gene_id: gene id in the .gff file.
                 type: int or str

    Returns: gene sequence in uppercase.
             type: str

    Dependencies: pandas as pd, get_genome_annotation(gff_path), get_gene_location(gff_path, gene_id), get_genome_sequence(fna_path), reverse_complement(seq)
    """
    contig, start, end, strand = get_gene_location(gff_path, gene_id)
    genome_seq = get_genome_sequence(fna_path)
    contig_seq = genome_seq[contig]
    if strand == '+':
        gene_seq = contig_seq[start - 1 : end]
    elif strand == '-':
        gene_seq = reverse_complement(contig_seq[start - 1 : end])
    return gene_seq

def get_snps(tsv_path, values=True):
    """Gets snps from a .tsv file.

    Args:
        tsv_path: path to the .tsv file with snps; only one letter for alternative allele is allowed.
                  type: str
        values: whether or not to read snp values,
                default is values=True.
                type: bool

    Returns: a table containing snp position, reference and alternative alleles, and snp values across a population.
             type: pd.dataframe

    Dependencies: pandas as pd
    """
    if values == False:
        snps = pd.read_csv(tsv_path, sep='\t', header=None, usecols=[1, 2, 3])
    else:
        snps = pd.read_csv(tsv_path, sep='\t', header=None)
    return snps