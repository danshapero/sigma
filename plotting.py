"""
Pseudocolor plots of solution to Poisson equation.

First command line argument is the name of the mesh, as stored in the
meshes/ folder, with no file extensions. For example, to use the mesh
described in the files
	meshes/circle.1.node
	meshes/circle.1.ele
the first argument would just be
	circle.1
	
Second command line argument is the file containing the solution, e.g.
	u.nc.
"""
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys
import math
import random
import numpy as np
from numpy import dtype
from scipy.io.netcdf import NetCDFFile as dataset


# Load the mesh files
nodes = open("meshes/"+sys.argv[1]+".node", "r")
n = int(nodes.readline().split()[0])
x = np.zeros(n)
y = np.zeros(n)
for i in range(n):
	line = nodes.readline().split()
	x[i] = float(line[1])
	y[i] = float(line[2])

elements = open("meshes/"+sys.argv[1]+".ele", "r")
m = int(elements.readline().split()[0])
ele = np.zeros((m,3))
for i in range(m):
	line = elements.readline().split()
	ele[i][:] = map(int,line[1:])
	for j in range(3):
		ele[i,j] = ele[i,j]-1

# Load the solution
ncfile = dataset(sys.argv[2], "r")
u = ncfile.variables['u'][:]
ncfile.close()

# Construct the triangulation
triang = tri.Triangulation(x,y,ele)

# Plot stuff
plt.figure()
plt.gca().set_aspect('equal')
plt.tricontourf(x, y, ele, u, 36, shading='faceted')
plt.colorbar()
plt.title('Solution of Lu = f')
plt.xlabel('x')
plt.ylabel('y')
plt.savefig(sys.argv[3]+'.png', dpi=100)
plt.clf()
#plt.show()
