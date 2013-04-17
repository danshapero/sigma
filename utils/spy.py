import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.sparse import *
from scipy import *

ia_file = open(sys.argv[1],"r")
ja_file = open(sys.argv[2],"r")
a_file = open(sys.argv[3],"r")

ia = array( map(int, ia_file.readlines() ) )
ja = array( map(int, ja_file.readlines() ) )
a = array( map(float, a_file.readlines() ) )

n = len(ia)-1
nnz = len(ja)

for i in range(n+1):
    ia[i] = ia[i]-1

for i in range(nnz):
    ja[i] = ja[i]-1

ia_file.close()
ja_file.close()
a_file.close()

A = csr_matrix( (a,ja,ia) )

plt.figure()
plt.spy(A,marker='.',markersize=0.5)
plt.savefig(sys.argv[4],dpi=100)
plt.clf()
