#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <graphs.h>

int main(int argc, char *argv[]) {

    ////////////////////////////////////////////////////////////////////////
    // Tests for C wrapper to graph data types                            //
    ////////////////////////////////////////////////////////////////////////

    // Create a graph on the Fortran side and get a C pointer to it
    void *cgp;

    get_graph(&cgp,0);
    graph_init(&cgp,100,100);

    int i = 0, j = 1;
    int con;
    con = -1;

    // Check whether two unconnected nodes actually are connected
    connected(&cgp,i,j,&con);
    if (con!=0) {
        printf("Nodes should not be connected for an empty graph\n");
    }

    // Add an edge and check to see that the nodes are connected
    add_edge(&cgp,i,j);
    connected(&cgp,i,j,&con);
    if (con!=1) {
        printf("Inserting edge unsuccessful\n");
    }

    // Delete an edge and check to see that the nodes are not connected
    delete_edge(&cgp,i,j);
    connected(&cgp,i,j,&con);
    if (con!=0) {
        printf("Deleting edge unsuccessful\n");
    }

    // Fill out the graph some more
    for (i=0; i<99; i++) {
        add_edge(&cgp,i,i+1);
        add_edge(&cgp,i+1,i);
    }

    // Find the neighbors of a given node and check that they are correct
    int *nbrs;
    nbrs = (int *)malloc( 4 * sizeof(int) );

    graph_neighbors(&cgp,50,nbrs,4);
    if ( !((nbrs[0]==51 && nbrs[1]==49) || (nbrs[0]==49 && nbrs[1]==51)) ) {
        printf("Getting neighbors of a node unsuccessful\n");
    }

    // Try compressing the graph to a more efficient format
    convert(&cgp,1);

    graph_neighbors(&cgp,50,nbrs,4);
    if ( !((nbrs[0]==51 && nbrs[1]==49) || (nbrs[0]==49 && nbrs[1]==51)) ) {
        printf("Getting neighbors of a node unsuccessful\n");
    }


    ////////////////////////////////////////////////////////////////////////
    // Tests for C wrapper to matrix data types                           //
    ////////////////////////////////////////////////////////////////////////
    void *cmp;

    get_matrix(&cmp,0);
    sparse_matrix_init(&cmp,100,100);




    return 0;

}
