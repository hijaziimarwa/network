import sys
import os



seed_min=int(sys.argv[1])
seed_incr=int(sys.argv[2])
seed_nsim=int(sys.argv[3])

fibers_min=int(sys.argv[4])
fibers_incr=int(sys.argv[5])
fibers_nsim=int(sys.argv[6])

for j in range(fibers_nsim):
    jj = fibers_min + fibers_incr*j
    for i in range(seed_nsim):
        ii = seed_min + seed_incr*i
        os.system('./change_seed.sh %d %d'%(ii, jj))
