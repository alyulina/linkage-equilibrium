import itertools as it
import pandas as pd
import numpy as np
import datetime


def gff_to_df(gff_path, save_path=None):
    """
    Parse GFF file downloaded from uhgg into dataframe
    GFF file contains predicted gene annotation
    E.g.:
        http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/uhgg_catalogue/MGYG-HGUT-000/MGYG-HGUT-00001/genome/MGYG-HGUT-00001.gff
    :param gff_path:
    :param save_path: If provided, save the pased file into a csv formatted file for future use
    :return:
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
    df = pd.DataFrame(list(zip(contigs, sources, types, starts, ends, scores, strands, frames, attributes)), columns=['Contig', 'Source', 'Type', 'Start', 'End', 'Score', 'Strand', 'Frame', 'Attribute'])
    df['Gene ID'] = np.arange(df.shape[0])  # We assign gene id by order in gff file. This is an arbitrary choice
    if save_path is not None:
        df.to_csv(save_path)
    return df


def process_SNVs_table(snvs_path, gene_df, save_path):
    """
    Split the SNV table into smaller ones by gene
    :param snvs_path: path to the big SNV table downloaded from uhgg and decompressed into .tsv
    :param gene_df: parsed by gff_to_df
    :param save_path: Path to folder holding all the individual gene files
    :return:
    """
    linecount = 0
    contig_dfs = {contig: df for contig, df in gene_df.groupby('Contig')}  # faster than slicing df every time changing contig
    with open(snvs_path, 'r') as f:
        head = f.readline()  # will be useful when processing individual genes
        for line in f:
            items = line[:-2].split('\t')
            contig = items[0]
            contig_df = contig_dfs[contig]
            pos = int(items[1])
            row = if_gene(contig_df, pos)  # sites found in multiple genes will be removed!
            if row is not None:
                # Important: here we set how the gene files should be named
                file_name = save_path + "{}-{}.tsv".format(contig, row['Gene ID'].squeeze())
                with open(file_name, 'a') as g:
                    g.write(line)
            if linecount % 100000 == 0:
                now = datetime.now()
                current_time = now.strftime("%H:%M:%S")
                print("Finished {} at {}".format(linecount, current_time))
            linecount += 1


def if_gene(df, pos):
    """
    Return the row containing information of the gene of a certain position
    it's possible for a position to be involved in multiple genes; we ignore those cases
    :param df: Gene df as parsed by gff_to_df
    :param pos: Position along the genome
    :return: If unique gene is found, return row; else None
    """
    mask = (df['Start'] <= pos) & (df['End'] >= pos)
    passed = np.sum(mask)
    if passed==1:
        return df[mask]
    else:
        return None

def get_non_redundant_genome_mask(mgnify_accession, genomes):
    """
    According to metadata, the total # of genomes can be a few times larger than the # of non-redundant genome
    Therefore, an additional step to filter genomes might be necessary
    :param mgnify_accession:
    :param genomes: List of genome names as in the header of SNV tables
    :return: a boolean mask for unique genomes
    """
    # TODO: currently doing nothing, but need to figure out how to deduplicate from metadata
    return np.ones(len(genomes)).astype(bool)


def process_SNV_table_single_gene(species_accession, header, filepath, gene_id, output=None):
    """
    Compute genotype counts for all pairs of sites in a single gene
    :param species_accession: mgnify ID
    :param header: the header for SNV table containing genome names etc; currently not in individual gene files
    :param filepath: file path to the gene SNV table
    :param gene_id: ID of the gene, assigned by us in the format "contig_id-gene_index"; useful to filter ncRNA or other genes later
    :param output: If provided, pairwise results will be written to the file as a line (n11, n10, n01, n00, ell) separated by whitespace
    :return: List of pairwise results, if output is not provided
    """

    dat = pd.read_csv(filepath, delimiter='\t', header=None)
    dat.columns = header.split('\t')
    dat.sort_values('Pos', inplace=True)

    # filtering sites with more than two alleles
    # TODO: this can be improved in the future to take care non-biallelic sites
    biallelic_dat = dat.groupby('Pos').filter(lambda x: x.shape[0] == 1)

    # converting to numpy array for faster arithmetics
    snvs = biallelic_dat.iloc[:, 4:].astype(int).to_numpy()
    genome_mask = get_non_redundant_genome_mask(species_accession, biallelic_dat.columns[4:])
    snvs = snvs[:, genome_mask]

    zeros = (snvs == 0).astype(int)
    ones = (snvs == 1).astype(int)
    n00s = np.dot(zeros, zeros.T)  # using dot product to count the number of genomes with 00 genotype
    n10s = np.dot(ones, zeros.T)
    n11s = np.dot(ones, ones.T)
    # TODO: in the future, we can also calculate pairs of syn/non-syn in a similar fashion here
    if output is not None:
        with open(output, 'a') as f:
            for i, j in it.combinations(range(snvs.shape[0]), 2):
                ell = np.abs(biallelic_dat.iloc[i, 1] - biallelic_dat.iloc[j, 1])
                n00 = str(n00s[i, j])
                n10 = str(n10s[i, j])
                n01 = str(n10s[j, i])
                n11 = str(n11s[i, j])
                res = (n11, n10, n01, n00, str(ell), gene_id)
                f.write(' '.join(res)+'\n')
    else:
        all_pairs = []
        for i, j in it.combinations(range(snvs.shape[0]), 2):
            ell = np.abs(biallelic_dat.iloc[i, 1] - biallelic_dat.iloc[j, 1])
            n00 = n00s[i, j]
            n10 = n10s[i, j]
            n01 = n10s[j, i]
            n11 = n11s[i, j]
            res = (n11, n10, n01, n00, ell, gene_id)
            all_pairs.append(res)
        return all_pairs
