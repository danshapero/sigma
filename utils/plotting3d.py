from mayavi import mlab
from tvtk.api import tvtk
import numpy as np
from scipy.io.netcdf import NetCDFFile as dataset
import sys

# Load the mesh files
nodes = open(sys.argv[1]+".node","r")
n = int(nodes.readline().split()[0])
x = np.zeros((n,3))
for i in range(n):
    x[i,:] = map(float,nodes.readline().split()[1:])

nodes.close()

elements = open(sys.argv[1]+".ele","r")
m = int(elements.readline().split()[0])
ele = np.zeros((m,4))
for i in range(m):
    ele[i,:] = map(int,elements.readline().split()[1:])
    for j in range(4):
        ele[i,j] = ele[i,j]-1

elements.close()

# Load the solution
ncfile = dataset(sys.argv[2],"r")
u = ncfile.variables['u'][:]
ncfile.close()

# Set up the VTK unstructured grid object
tet_type = tvtk.Tetra().cell_type
grid = tvtk.UnstructuredGrid(points=x)
grid.set_cells(tet_type,ele)

# Add the scalar field data to the grid
grid.point_data.scalars = u
grid.point_data.scalars.name = "u"

# Some nonsense with the Mayavi pipeline
ds = mlab.pipeline.add_dataset(grid)
iso = mlab.pipeline.iso_surface(ds)
iso.actor.property.opacity = 0.5
iso.contour.number_of_contours = 10
mlab.show()
