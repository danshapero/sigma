import sys
import numpy as np

def read_poly_file(filename):
    # Read in the .poly file specified in the first command line argument
    fid = open(filename,"r")
    n = int(fid.readline().split()[0])

    x = np.zeros((n,2))
    for i in range(n):
        x[i,:] = map(float,fid.readline().split()[1:3])

    m = int(fid.readline().split()[0])
    lines = np.zeros((m,2))
    for i in range(m):
        lines[i,:] = map(int,fid.readline().split()[1:3])
    lines = lines.astype(int)

    p = int(fid.readline().split()[0])
    holes = np.zeros((p,2))
    for i in range(p):
        holes[i,:] = map(float,fid.readline().split()[1:3])
    fid.close()

    return (x,lines,holes)


def find_loops(lines):
    m = len(lines)
    n = int(np.max(lines))

    # Generate the requisite line loops needed by gmsh
    next_node = range(n)
    for i in range(m):
        next_node[ lines[i,0]-1 ] = int(lines[i,1]-1)

    loop = []
    current_loop = -1
    nodes_not_entered_yet = range(n)
    while nodes_not_entered_yet:
        current_loop = current_loop+1
        loop.append( [] )
        head = nodes_not_entered_yet.pop(0)
        loop[current_loop].append( head )
        tail = head
        while ( next_node[tail] != head ):
            tail = next_node[tail]
            loop[current_loop].append( tail )
            nodes_not_entered_yet.pop( nodes_not_entered_yet.index(tail) )
        loop[current_loop].append( head )

    return loop


if __name__ == "__main__":
    (x,lines,holes) = read_poly_file(sys.argv[1])

    n = len(x)
    m = len(lines)
    p = len(holes)

    loop = find_loops(lines)
    l = len(loop)

    fid = open(sys.argv[2],"w")
    fid.write("lc = 1e-2\n")
    for i in range(n):
        fid.write("Point({0}) = {{ {1}, {2}, lc }} ;\n".format(i+1,x[i,0],x[i,1]))

    fid.write("\n")
    sm = 1
    for i in range(l):
        q = len(loop[i])
        for j in range(q-1):
            fid.write("Line({0}) = {{ {1},{2} }} ;\n".format(j+sm,loop[i][j]+1,loop[i][j+1]+1))
        sm = sm+q-1

    fid.write("\n")
    sm = 1
    for i in range(l):
        q = len(loop[i])
        fid.write("Line Loop({0}) = {{ ".format(m+i+1))
        for j in range(q-2):
            fid.write("{0}, ".format(j+sm))
        fid.write("{0} }} ;\n".format(sm+q-2))
        sm = sm+q-1

    fid.write("\n")
    fid.write("Plane Surface({0}) = {{ ".format(m+l+1))
    for i in range(l-1):
        fid.write("{0}, ".format(m+i+1))
    fid.write("{0} }} ;\n".format(m+l))
