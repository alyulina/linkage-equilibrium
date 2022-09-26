import os
from typing import List

import genome_utils
from optparse import OptionParser

#ref_genome_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/MGYG-HGUT-02506.fna"
#genome_annotation_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/MGYG-HGUT-02506.gff"
#snps_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test"

#ann_snps_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test"
#err_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/MGYG-HGUT-02506_snp_annotation.err"

parser = OptionParser()

parser.add_option("-s",
                  "--species",
                  type="string",
                  help="species uhgg accession, e. g. MGYG-HGUT-00001",
                  dest="species_accession")

(options, args) = parser.parse_args()
species_accession = options.species_accession

ref_genome_path = '/home/groups/bhgood/uhgg/reference_genomes/' + species_accession + '.fna'
genome_annotation_path = '/home/groups/bhgood/uhgg/genes/' + species_accession + '.gff'
snps_path = '/scratch/groups/bhgood/uhgg/snv_tables/' + species_accession # path to folder

ann_snps_path = '/scratch/groups/bhgood/uhgg/snv_tables_annotated/' + species_accession # path to folder
err_path = '/home/users/alyulina/recombination/uhgg/' + species_accession + '_snp_annotation.err'

err_out = open(err_path, 'a')

if os.path.exists(ann_snps_path) == False:
    os.makedirs(ann_snps_path)

ref_genome = genome_utils.get_genome_sequence(ref_genome_path)
genome_annotation = genome_utils.get_genome_annotation(genome_annotation_path)

for i in range(len(genome_annotation)): # for each gene

    snp_path = snps_path + '/' + genome_annotation.iloc[i, 0] + '-' + str(genome_annotation.iloc[i, 9]) + '.tsv'
    if os.path.exists(snp_path) == False: # skipping the gene if it does not have snps
        continue

    if genome_annotation.iloc[i, 2] != 'CDS': # skipping noncoding genes
        continue

    gene_id = genome_annotation.iloc[i, 9]
    contig, start, end, strand = genome_utils.get_gene_location(genome_annotation_path, gene_id) # numbering starts w/ 1
    gene_seq = genome_utils.get_gene_sequence(ref_genome_path, genome_annotation_path, gene_id)
    snps = genome_utils.get_snps(snp_path)
    n_snps = len(snps)

    err_out.write('gene id: ' + str(gene_id) + '\n')

    site_degeneracies = [0 for x in range(n_snps)] # default for ambiguous snps
    snp_types = ['-' for x in range(n_snps)] # default for ambiguous snps

    snp_type_counts = [0, 0, 0]  # s, n, m
    n_not_snps = 0
    for j in range(n_snps):

        # checking that there are no ambiguous nucleotides
        if (snps.iloc[j, 2] not in ['A', 'T', 'G', 'C']) or (snps.iloc[j, 3] not in ['A', 'T', 'G', 'C']):
            err_out.write('ambiguous snp in line ' + str(j) + '\n')
            continue

        snp_pos_genome = snps.iloc[j, 1] - 1

        if strand == '+':
            snp_pos_gene = snp_pos_genome - (start - 1)
            ref_allele = snps.iloc[j, 2]
            alt_allele = snps.iloc[j, 3]
        elif strand == '-':
            snp_pos_gene = (end - 1) - snp_pos_genome
            ref_allele = genome_utils.reverse_complement(snps.iloc[j, 2])
            alt_allele = genome_utils.reverse_complement(snps.iloc[j, 3])
        else:
            continue

        if ref_allele != gene_seq[snp_pos_gene]:
            err_out.write('snps do not map in line ' + str(j) + '\n')

        codon_pos = snp_pos_gene // 3 # codon loc within gene
        snp_pos_codon = snp_pos_gene % 3 # snp loc within codon
        ref_codon = gene_seq[codon_pos * 3 : (codon_pos + 1) * 3]
        alt_codon = ref_codon[:snp_pos_codon] + alt_allele + ref_codon[(snp_pos_codon + 1):]
        #print(ref_codon, alt_codon)

        #if ref_codon not in genome_utils.codon_table: # usually because reference has an N; skipping it
        #    err_out.write('skipped an ambiguous ref codon in line ' + str(j) + '\n')
        #    continue

        #if alt_codon not in genome_utils.codon_table: # usually because reference has an N; skipping it
        #    err_out.write('skipped an ambiguous alt codon in line ' + str(j) + '\n')
        #    continue

        site_degeneracies[j] = genome_utils.codon_degeneracy_table[ref_codon][snp_pos_codon]

        if ref_codon == alt_codon:
            snp_type = '.'
            n_not_snps += 1
        elif genome_utils.codon_table[ref_codon] == genome_utils.codon_table[alt_codon]:
            snp_type = 's' # synonymous
            snp_type_counts[0] += 1
        elif genome_utils.codon_table[ref_codon] == '!' or genome_utils.codon_table[alt_codon] == '!':
            # premature gain or loss of stop codon
            snp_type = 'n' # nonsense
            snp_type_counts[1] += 1
        else:
            snp_type = 'm'  # missense
            snp_type_counts[2] += 1

        snp_types[j] = snp_type

    if n_not_snps != 0:
        err_out.write(str(n_not_snps) + ' out of ' + str(len(snps)) + ' snps are not snps\n')

    snps.insert(4, 'snp type', snp_types)
    snps.insert(5, 'site degeneracy', site_degeneracies)
    snps.to_csv(ann_snps_path + '/' + genome_annotation.iloc[i, 0] + '-' + str(gene_id) + '-annotated.tsv', sep='\t', header=False, index=False)
    #err_out.write('synonymous: ' + str(snp_type_counts[0]) + ', nonsynonumous: ' + str(snp_type_counts[1]) + ', missense: ' + str(snp_type_counts[2]) + '\n')
    err_out.write('done!\n\n')

err_out.close()
