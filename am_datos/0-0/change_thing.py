import numpy as np

dat = np.loadtxt('magnetic.dat')
dat[:,-1] = 0/dat[:,-1]
np.savetxt('out.txt', dat, fmt='%.2f')
