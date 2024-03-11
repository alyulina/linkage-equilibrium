import numpy as np
import gzip
import os
import argparse

import parameters

parser = argparse.ArgumentParser()

# Add the arguments
parser.add_argument("-p", "--parameters", type=str, 
                    help="name of the parameters regime from parameters.py", 
                    dest="regime")
parser.add_argument("--path", type=str, 
                    help="root path to the data files", 
                    dest="dat_path")
parser.add_argument("--debug", action='store_true', 
                    help="Enable debug mode")
parser.add_argument("--save-path", type=str, 
                    help="root path to where npy files will be saved", 
                    dest="savepath")

# Parse the arguments
args = parser.parse_args()

# Access the arguments
regime = args.regime
dat_path = args.dat_path
save_path = args.savepath
debug = args.debug

params = parameters.params[regime]  # list of tuples of parameters
rs = np.array([params[idx][6] for idx in range(len(params))])

for param_idx in range(len(rs)):
    filename = os.path.join(dat_path, 'output_%s_%d.txt.gz' % (regime,param_idx))

    file = gzip.GzipFile(filename,"r")
    f11s = []
    f10s = []
    f01s = []

    for line in file:

        if line.startswith(b'//'):
            continue
        items = line.split()

        f11 = float(items[0])
        f10 = float(items[1])
        f01 = float(items[2])

        f11s.append(f11)
        f10s.append(f10)
        f01s.append(f01)

    file.close()

    f11s = np.array(f11s)
    f10s = np.array(f10s)
    f01s = np.array(f01s)
    f00s = 1-f11s-f10s-f01s
    fAs = f11s+f10s
    fBs = f11s+f01s
    
    all_fs = np.stack((f11s, f10s, f01s), axis=-1)
    np.save(os.path.join(save_path, 'fs_%s_%d'%(regime,param_idx)), all_fs)
    if debug:
        print("debugging: exiting after one rep")
        break
