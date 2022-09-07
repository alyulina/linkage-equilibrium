import os
import genome_utils

#ref_genome_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/MGYG-HGUT-02492.fna"
#genome_annotation_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/MGYG-HGUT-02492.gff"
#snps_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test"

#ann_snps_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test"
#err_path = "/Users/alyulina/Projects/Linkage equilibrium/uhgg/test/MGYG-HGUT-02492_snp_annotation.err"

ref_genome_path = '/home/groups/bhgood/uhgg/reference_genomes/MGYG-HGUT-02492.fna'
genome_annotation_path = '/home/groups/bhgood/uhgg/genes/MGYG-HGUT-02492.gff'
snps_path = '/home/groups/bhgood/uhgg/snv_tables/MGYG-HGUT-02492' # path to folder

ann_snps_path = '/home/groups/bhgood/uhgg/snv_tables_annotated' # path to folder
err_path = '/home/users/alyulina/recombination/uhgg/MGYG-HGUT-02492_snp_annotation.err'
err_out = open(err_path, 'a')

ref_genome = genome_utils.get_genome_sequence(ref_genome_path)
genome_annotation = genome_utils.get_genome_annotation(genome_annotation_path)

for i in range(len(genome_annotation)):

    snp_path = snps_path + '/' + genome_annotation.iloc[i, 0] + '-' + str(genome_annotation.iloc[i, 9]) + '.tsv'
    if os.path.exists(snp_path) == False:
        continue

    if genome_annotation.iloc[i, 2] != 'CDS':
        continue

    gene_id = genome_annotation.iloc[i, 9]
    start, end, strand = genome_utils.get_gene(genome_annotation_path, gene_id) # numbering starts w/ 1
    gene_seq = genome_utils.get_gene_sequence(ref_genome_path, genome_annotation_path, gene_id)
    snps = genome_utils.get_snps(snp_path)

    err_out.write('gene id: ' + str(gene_id) + '\n')

    site_degeneracies = []
    snp_types = []

    n_not_snps = 0
    for j in range(len(snps)):
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

        codon_pos = snp_pos_gene // 3
        snp_pos_codon = snp_pos_gene % 3
        ref_codon = gene_seq[codon_pos * 3 : (codon_pos + 1) * 3]
        alt_codon = ref_codon[:snp_pos_codon] + alt_allele + ref_codon[(snp_pos_codon + 1):]
        # print(ref_codon, alt_codon)

        site_degeneracies.append(genome_utils.codon_degeneracy_table[ref_codon][snp_pos_codon])

        if ref_codon == alt_codon:
            snp_type = '.'
            n_not_snps += 1
        elif genome_utils.codon_table[ref_codon] == genome_utils.codon_table[alt_codon]:
            snp_type = 's' # synonymous
        elif genome_utils.codon_table[ref_codon] == '!' or genome_utils.codon_table[alt_codon] == '!':
            # premature gain or loss of stop codon
            snp_type = 'n' # nonsense
        else:
            snp_type = 'm'  # missense

        snp_types.append(snp_type)

        if n_not_snps != 0:
            err_out.write(str(n_not_snps) + ' out of ' + str(len(snps)) + ' snps are not snps\n')

    snps.insert(4, 'snp type', snp_types)
    snps.insert(5, 'site degeneracy', site_degeneracies)
    snps.to_csv(ann_snps_path + '/' + genome_annotation.iloc[i, 0] + '-' + str(gene_id) + '-annotated.tsv',
                sep='\t', header=False, index=False)
    err_out.write('done!\n\n')

err_out.close()




