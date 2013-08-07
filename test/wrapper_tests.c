#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[]) {

    void *cgp;

    get_graph(&cgp,0);
    graph_init(&cgp,100,100);
    yea_bitches(&cgp);

    return 0;

}
