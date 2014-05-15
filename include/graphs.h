// C header file for the routines defined in wrapper.f03

void get_graph(void *cgp, int storage_format);
void graph_init(void *cgp, int n, int m);
void graph_neighbors(void *cgp, int *nbrs, int i, int d);
void connected(void *cgp, int i, int j, int *con);
void add_edge(void *cgp, int i, int j);
void delete_edge(void *cgp, int i, int j);
