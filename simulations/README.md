This folder contains [Ben's code from his 2022 Genetics paper](https://github.com/bgoodlab/rare_ld) + Zhiru's code for the new statistic. Here's how to use it.  
  
First, compile the C++ code: `g++ -o simulate_twolocus twolocus.cpp`.  Do not change the executable name.  Use g++ or another compiler installed. On the cluster, type `g++ -o simulate_twolocus -std=c++11 twolocus.cpp`.   

Create the outpute folder: `mkdir output`. The name can be changed in the bash script. If running from /home on the cluster, make sure that the output is stored in /scratch.    
  
Then run the bash script that will run simulations: `bash run_simulation.sh r`, where r specifies the parameter regime from parameters.py. Or, on the cluster, `sbatch srun_simulations.sh r`. The python code needs python 2, so be sure to load it if running on the cluster or change python to python2 in the bash script if running locally.  
  
To plot the linkage equilibrium statistic locally, use `python2 plot_le_recombination.py -p r`. If running on the cluster, use `sbatch srun_plot_le_recombination.sh r`. This script uses ld_theory.py, which requires scipy version 0.17.0 or higher. If running locally, switch to the conda environment named ben: `conda activate ben`. If running on the cluster, make sure to locally install all dependencies.  
  
Note: make sure to change this code to the version w/ recurrent muts + check that compute_le is up to date
