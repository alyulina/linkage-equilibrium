import numpy
import sys

params = {}

num_runs = int(1e9)
N = 1e06
dt = 100
s1 = 0
s2 = 0
eps = 0
rs = numpy.hstack([[0],numpy.logspace(-8,-1,21)])

params['test'] = [(100000,dt,N,s1,s2,eps,r) for r in rs]

# type == 'r'
params['r'] = [(num_runs,dt,N,s1,s2,eps,r) for r in rs]
params['r_small'] = [(int(1e6),dt,N,s1,s2,eps,r) for r in rs]

mu = 1e-10
params['r_mu=1e-10'] = [(num_runs,dt,N,s1,s2,eps,r,mu) for r in rs]

mu = 1e-9
params['r_mu=1e-9_small'] = [(int(1e6),dt,N,s1,s2,eps,r,mu) for r in rs]
params['r_mu=1e-9'] = [(num_runs,dt,N,s1,s2,eps,r,mu) for r in rs]

mu = 1e-8
params['r_mu=1e-8_small'] = [(int(1e6),dt,N,s1,s2,eps,r,mu) for r in rs]
params['r_mu=1e-8'] = [(num_runs,dt,N,s1,s2,eps,r,mu) for r in rs]


s = 1e-5
params['r_selA=1e-5'] = [(num_runs,dt,N,s,0,0,r) for r in rs]
params['r_selAB=1e-5'] = [(num_runs,dt,N,s/2,s/2,0,r) for r in rs]
params['r_eps=1e-5'] = [(num_runs,dt,N,0,0,s,r) for r in rs]
params['r_negeps=1e-5'] = [(num_runs,dt,N,s,s,-s,r) for r in rs]
params['r_negeps_strong=1e-5'] = [(num_runs,dt,N,s,s,-2*s,r) for r in rs]

s = 1e-4
params['r_selA=1e-4'] = [(num_runs,dt,N,s,0,0,r) for r in rs]
params['r_selAB=1e-4'] = [(num_runs,dt,N,s/2,s/2,0,r) for r in rs]
params['r_eps=1e-4'] = [(num_runs,dt,N,0,0,s,r) for r in rs]
params['r_negeps=1e-4'] = [(num_runs,dt,N,s,s,-s,r) for r in rs]
params['r_negeps_strong=1e-4'] = [(num_runs,dt,N,s,s,-2*s,r) for r in rs]

s = 1e-3
params['r_selA=1e-3'] = [(num_runs,dt,N,s,0,0,r) for r in rs]
params['r_selAB=1e-3'] = [(num_runs,dt,N,s/2,s/2,0,r) for r in rs]
params['r_eps=1e-3'] = [(num_runs,dt,N,0,0,s,r) for r in rs]
params['r_negeps=1e-3'] = [(num_runs,dt,N,s,s,-s,r) for r in rs]
params['r_negeps_strong=1e-3'] = [(num_runs,dt,N,s,s,-2*s,r) for r in rs]

s = 1e-2
params['r_selA=1e-2'] = [(num_runs,dt,N,s,0,0,r) for r in rs]
params['r_selAB=1e-2'] = [(num_runs,dt,N,s/2,s/2,0,r) for r in rs]
params['r_eps=1e-2'] = [(num_runs,dt,N,0,0,s,r) for r in rs]
params['r_negeps=1e-2'] = [(num_runs,dt,N,s,s,-s,r) for r in rs]
params['r_negeps_strong=1e-2'] = [(num_runs,dt,N,s,s,-2*s,r) for r in rs]


ss = numpy.logspace(-5,-2,13)
params['selA'] = [(num_runs,dt,N,s,0,0,0) for s in ss]
params['selB'] = [(num_runs,dt,N,0,s,0,0) for s in ss]
params['selAB'] = [(num_runs,dt,N,s/2,s/2,0,0) for s in ss]
params['eps'] = [(num_runs,dt,N,0,0,s,0) for s in ss]
params['negeps'] = [(num_runs,dt,N,s,s,-s,0) for s in ss]

# params['selAB_v_small_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[0]) for s in ss] # selAB + very small r
# params['selAB_small_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[1]) for s in ss] # selAB + small r
# params['selAB_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[2]) for s in ss] # selAB + r
# params['selAB_large_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[3]) for s in ss] # selAB + large r
# params['selAB_v_large_r'] = [(num_runs,dt,N,s/2,s/2,0,numpy.logspace(-8,-1,5)[4]) for s in ss] # selAB + very large r


if __name__=='__main__':

    if sys.argv[1]=='idxs':
        type = sys.argv[2].strip()
        for r_idx in xrange(0,len(params[type])):
            print r_idx
    elif sys.argv[1]=='get_params':
        idx = long(sys.argv[3])
        type = sys.argv[2].strip()
        print " ".join([str(item) for item in params[type][idx]])
