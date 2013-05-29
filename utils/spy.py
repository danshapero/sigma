import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.sparse import *
from scipy import *

fid = open(sys.argv[1],"r")

line = fid.readline().split()
(nrow,ncol,nnz) = map(int,line)

rows = np.zeros(nnz,int)
cols = np.zeros(nnz,int)
a = np.zeros(nnz)

for i in range(nnz):
    line = fid.readline().split()
    if not line: break
    rows[i] = int(line[0])
    cols[i] = int(line[1])
    a[i] = float(line[2])

fid.close()

for i in range(nnz):
    rows[i] = rows[i]-1
    cols[i] = cols[i]-1

A = coo_matrix( (a,(rows,cols)) )

plt.figure()
plt.spy(A,marker='.',markersize=0.5)
plt.savefig(sys.argv[2],dpi=100)
plt.clf()
