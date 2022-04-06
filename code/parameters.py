import numpy
import sys

params = {}

num_runs = 10000000
# num_runs = 100000000
N = 1e05
dt = 100
s1 = 0
s2 = 0
eps = 0
rs = numpy.hstack([[0],numpy.logspace(-8,-1,21)])

# type == 'r'
params['r'] = [(num_runs,dt,N,s1,s2,eps,r) for r in rs]
params['r_s=1e-5'] = [(num_runs,dt,N,0.5e-5,0.5e-5,eps,r) for r in rs]
params['r_s=1e-4'] = [(num_runs,dt,N,0.5e-4,0.5e-4,eps,r) for r in rs]
params['r_s=1e-3'] = [(num_runs,dt,N,0.5e-3,0.5e-3,eps,r) for r in rs]
params['r_s=1e-2'] = [(num_runs,dt,N,0.5e-2,0.5e-2,eps,r) for r in rs]
params['r_eps=1e-5'] = [(num_runs,dt,N,0,0,1e-5,r) for r in rs]
params['r_eps=1e-4'] = [(num_runs,dt,N,0,0,1e-4,r) for r in rs]
params['r_eps=1e-3'] = [(num_runs,dt,N,0,0,1e-3,r) for r in rs]
params['r_eps=1e-2'] = [(num_runs,dt,N,0,0,1e-2,r) for r in rs]

rs = numpy.hstack([[0],numpy.logspace(-9,-4,21)])
params['small_r'] = [(num_runs,dt,N,s1,s2,eps,r) for r in rs]

ss = numpy.logspace(-5,-2,13)
params['selA'] = [(num_runs,dt,N,s,0,0,0) for s in ss]
params['selB'] = [(num_runs,dt,N,0,s,0,0) for s in ss]
params['selAB'] = [(num_runs,dt,N,s/2,s/2,0,0) for s in ss]
params['eps'] = [(num_runs,dt,N,0,0,s,0) for s in ss]
params['negeps'] = [(num_runs,dt,N,s,s,-s,0) for s in ss]

params['selAB_v_small_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[0]) for s in ss] # selAB + very small r
params['selAB_small_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[1]) for s in ss] # selAB + small r
params['selAB_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[2]) for s in ss] # selAB + r
params['selAB_large_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[3]) for s in ss] # selAB + large r
params['selAB_v_large_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[4]) for s in ss] # selAB + very large r

s1 = 1e-03/8
s2 = 1e-03/8
eps = 0 
rs = numpy.hstack([[0],numpy.logspace(-6,-2,17)])
params['r_selA'] = [(num_runs,dt,N,s1,s2,eps,r) for r in rs]


if __name__=='__main__':

    if sys.argv[1]=='idxs':
        type = sys.argv[2].strip()
        for r_idx in xrange(0,len(params[type])):
            print r_idx
    elif sys.argv[1]=='get_params':
        idx = long(sys.argv[3])
        type = sys.argv[2].strip()
        print " ".join([str(item) for item in params[type][idx]])
