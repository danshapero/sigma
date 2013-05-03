import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.sparse import *
from scipy import *

fid = open(sys.argv[1],"r")

line = fid.readline().split()
n = int(line[0])
nnz = int(line[1])

rows = np.zeros(nnz,int)
cols = np.zeros(nnz,int)
a = np.zeros(nnz)

i = 0
while True:
    line = fid.readline().split()
    if not line: break
    if len(line)==1:
        k = int(line[0])
    else:
        rows[i] = k
        cols[i] = int(line[0])
        a[i] = float(line[1])
        i = i+1

fid.close()

for i in range(nnz):
    rows[i] = rows[i]-1
    cols[i] = cols[i]-1

A = coo_matrix( (a,(rows,cols)) )

plt.figure()
plt.spy(A,marker='.',markersize=0.5)
plt.savefig(sys.argv[2],dpi=100)
plt.clf()
