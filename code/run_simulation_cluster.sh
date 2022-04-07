
#!/bin/bash
#
#This is a job array to run Ben's code 10 times for a set of parameters from parameters.py
#
#Give your job a name
#SBATCH --job-name=LE_r
#SBATCH --output=/scratch/users/alyulina/recombination/output/slurm-%A_%a.out
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
#SBATCH --array=1-10
#
#Do not export the local environment to the compute nodes
#SBATCH --export=NONE
unset SLURM_EXPORT_ENV
#
#For single-CPU jobs, make sure that they use a single thread
export OMP_NUM_THREADS=1
#
#Load the software module you intend to use
module load python/python/2.7.13
module load py-numpy/1.14.3_py27
module load py-scipy/1.1.0_py27 # need 0.17.0 or higher
#
export type=$1
for idx in $(python2 parameters.py idxs ${type}); do
    export params=$(python2 parameters.py get_params ${type} ${idx})
    echo ${params}
    ./simulate_twolocus ${params} | gzip -c > output/output_${type}_${idx}_%a.txt.gz
done