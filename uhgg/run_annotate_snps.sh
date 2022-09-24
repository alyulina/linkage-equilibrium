#!/bin/bash
#
#SBATCH --job-name=snps_ann
#SBATCH --output=/home/users/alyulina/recombination/uhgg/slurm-%j.out
#SBATCH --time=2-00:00:00
#SBATCH --mem=8G
#SBATCH --partition hns,dpetrov,normal
#SBATCH --mail-user=alyulina@stanford.edu
#SBATCH --mail-type=NONE
#SBATCH --no-requeue
#SBATCH --array=1-10
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
export OMP_NUM_THREADS=1
#
module load python/3.6.1
module load py-numpy/1.18.1_py36
module load viz
module load py-pandas/1.0.3_py36
#
echo $SLURM_ARRAY_TASK_ID
srun --cpu_bind=verbose python3 annotate_snps.py -s $(head -$SLURM_ARRAY_TASK_ID species_accessions.txt | tail -1)
