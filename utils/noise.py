import sys
import math
import random
import numpy as np
from numpy import dtype
from scipy.io.netcdf import NetCDFFile as dataset


nodes = open(sys.argv[1]+".node","r")
n = int(nodes.readline().split()[0])
nodes.close()

f = np.zeros(n)

sigma = 1.0
if len(sys.argv) == 4:
    sigma = float(sys.argv[3])
for i in range(n):
    f[i] = random.gauss(0,sigma)

ncfile = dataset(sys.argv[2],"w")
ncfile.createDimension("nodes",n)
data = ncfile.createVariable("u",dtype("float32").char,("nodes",))
data[:] = f[:]
ncfile.close()
