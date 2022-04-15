#!/bin/bash
#
#This is a job array to run Ben's code for a set of parameters from parameters.py
#
#Give your job a name
#SBATCH --job-name=LE
#SBATCH --output=/home/users/alyulina/recombination/output/slurm-%A_%a.out
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
#SBATCH --mail-type=ALL
#
#Do not restart the job if it fails
#SBATCH --no-requeue
#
#Submit a job array of N jobs (N is limited to 1000)
#Use the first and the last indices form $(python2 parameters.py idxs ${type})
#SBATCH --array=0-21
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
#
#Read type from parameters.py and print out parameters, then run simulation
type=$1
echo $(python2 parameters.py get_params ${type} $SLURM_ARRAY_TASK_ID)
srun --cpu_bind=verbose ./simulate_twolocus $(python2 parameters.py get_params ${type} $SLURM_ARRAY_TASK_ID) | gzip -c > /scratch/users/alyulina/recombination/output/output_${type}_$SLURM_ARRAY_TASK_ID.txt.gz
