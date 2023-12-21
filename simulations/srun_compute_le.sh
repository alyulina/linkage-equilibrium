#!/bin/bash
#
#Give your job a name
#SBATCH --job-name=redo_2022_figs
#SBATCH --output=/home/users/zhiru/linkage-equillibrium/code/batch_out/slurm-%j.out
#SBATCH --error=/home/users/zhiru/linkage-equillibrium/code/batch_out/slurm-%j.err
#
#Specify time limit; days-hours:minutes:seconds or hours:minutes:seconds
#SBATCH --time=0-12:00:00
#
#Specify memory in gigabytes
#SBATCH --mem=64G
#
#Specify account and partition
#SBATCH --partition normal
#
#Would you like to be notified when the job starts or is completed?
#SBATCH --mail-user=zhiru@stanford.edu
#SBATCH --mail-type=ALL
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
ml py-scipystack/1.0_py27
#
#
type=$1
srun --cpu_bind=verbose python2 compute_simulated_le_stats.py -p ${type} --path /scratch/groups/bhgood/LE_data/
#srun --cpu_bind=verbose python2 plot_recombination_figure.py
#srun --cpu_bind=verbose python2 plot_le_recombination.py -p ${type} --path $SCRATCH
