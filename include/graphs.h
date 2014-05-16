// C header file for the routines defined in wrapper.f03

void get_graph(void *cgp, int storage_format);
void graph_init(void *cgp, int n, int m);
void degree(void *cgp, int i, int *d);
void get_neighbors(void *cgp, int *nbrs, int i, int d);
void connected(void *cgp, int i, int j, int *con);
void add_edge(void *cgp, int i, int j);
void delete_edge(void *cgp, int i, int j);
void left_permute(void *cgp, int *p, int n);
void right_permute(void *cgp, int *p, int m);
