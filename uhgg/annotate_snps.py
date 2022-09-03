import genome_utils

ref_genome_path = '/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/MGYG-HGUT-02492.fna'
genome_annotation_path = '/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/MGYG-HGUT-02492.gff'
gene_path = '/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/GUT_GENOME143712_1-0.tsv'

ref_genome = genome_utils.get_genome_sequence(ref_genome_path)
genome_annotation = genome_utils.get_genome_annotation(genome_annotation_path)
# print(genome_annotation)
# get a list of genes from genome annotation, write out to std err if there's no file
# else get number of genes, snps paths, and gene ids

# have a for loop here
gene_id_i = 0 # same format as in genome_annotation
gene_seq_i = genome_utils.get_gene_sequence(ref_genome_path, genome_annotation_path, gene_id_i)

# if set(gene_seq_i) != {'C', 'A', 'T', 'G'}:
# print err in std err

snps_path_i = '/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/GUT_GENOME143712_1-' + str(gene_id_i) + '.tsv'

# add two empty columns

# go line by line:
# snp_id = ... (get first column value)
# determine 1) codon based on position (position // 3) and location within codon (position % 3)
# copy from Ben's code
codon_id = snp_id // 3
codon_position = snp_id % 3
codon = gene_seq[codon_id * 3 : (codon_id + 1) * 3]
genome_base = .. # - ref nuc us correct, write to std err if not; have a counter?
# - alt nuc differs from ref, write to std err if not; have a counter?


