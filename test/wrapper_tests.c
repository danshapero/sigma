#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <graphs.h>

int main(int argc, char *argv[]) {

    void *cgp;

    get_graph(&cgp,0);
    graph_init(&cgp,100,100);
    int i = 1, j = 2;
    int con;
    con = -1;

    connected(&cgp,i,j,&con);
    if (con!=0) {
        printf("Nodes should not be connected for an empty graph\n");
    }
    add_edge(&cgp,i,j);
    connected(&cgp,i,j,&con);
    if (con!=1) {
        printf("Inserting edge unsuccessful\n");
    }
    delete_edge(&cgp,i,j);
    connected(&cgp,i,j,&con);
    if (con!=0) {
        printf("Deleting edge unsuccessful\n");
    }

    return 0;

}
