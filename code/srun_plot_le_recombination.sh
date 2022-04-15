#!/bin/bash
#
#Give your job a name
#SBATCH --job-name=plot_LE
#SBATCH --output=/home/users/alyulina/recombination/output/slurm-%j.out
#
#Specify time limit; days-hours:minutes:seconds or hours:minutes:seconds
#SBATCH --time=2-00:00:00
#
#Specify memory in gigabytes
#SBATCH --mem=4G
#
#Specify account and partition
#SBATCH --partition hns,dpetrov,normal
#
#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=alyulina@stanford.edu
#SBATCH --mail-type=NONE
#
#Do not restart the job if it fails
#SBATCH --no-requeue
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#For single-CPU jobs, make sure that they use a single thread
export OMP_NUM_THREADS=1
#
#Load the software module you intend to use
module load python/2.7.13
module load py-numpy/1.14.3_py27
module load py-scipy/1.1.0_py27
module load viz
module load py-matplotlib/2.1.2_py27
#
#
type=$1
srun --cpu_bind=verbose python2 plot_le_recombination.py -p ${type}