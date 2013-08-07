#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {

    void *cgp;

    get_graph(&cgp,0);
    graph_init(&cgp,100,100);
    int i = 1, j = 2;
    int con;
    con = -1;

    connected(&cgp,i,j,&con);
    printf("%d\n",con);
    add_edge(&cgp,i,j);
    connected(&cgp,i,j,&con);
    printf("%d\n",con);
    delete_edge(&cgp,i,j);
    connected(&cgp,i,j,&con);
    printf("%d\n",con);

    return 0;

}
