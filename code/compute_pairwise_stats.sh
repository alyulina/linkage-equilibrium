#!/usr/bin/bash
#SBATCH --job-name=compute_half
#SBATCH --output=out/compute_half_%j.%a.out
#SBATCH --error=out/compute_half_%j.%a.err
#SBATCH --time=18:00:00
#SBATCH -p normal 
#SBATCH -c 1
#SBATCH --mem=64GB
#SBATCH --array=1-1

module load py-pandas/1.0.3_py36
module load py-sympy/1.1.1_py36
module load py-scipy/1.1.0_py36
export PYTHONPATH=$PYTHONPATH:"${HOME}/linkage-equillibrium/code/"

#while read p; do
#  accession=$p
#  filedir="${GROUP_HOME}/uhgg/site_pairs/${accession}/${accession}.npy"
#  python3 ../compute_rare_LE_microbiome.py LD --path $filedir --savepath $savedir
#done < "accession_list.txt"

# to use slurm array to compute all 10 species
#accession=`sed -n "${SLURM_ARRAY_TASK_ID}p" < accession_list.txt`
accession="MGYG-HGUT-02492"

filedir="${GROUP_HOME}/uhgg/site_pairs/${accession}/${accession}.npy"
savedir="${GROUP_HOME}/uhgg/pairwise_stats/"
#filedir="./erectale_snv_pairs.txt.gz"
#savedir="./test_pairwise_stats/"
#python3 ../compute_rare_LE_microbiome.py LE --path $filedir --savepath $savedir --pairtype 0
#python3 ../compute_rare_LE_microbiome_alternative.py LE --path $filedir --savepath $savedir --pairtype 0
python3 ../compute_half_rare_LE_microbiome.py --path $filedir --savepath $savedir --pairtype 0
