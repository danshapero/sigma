extern"C" {
typedef struct {
    void *p;
    int nrow, ncol, nnz, max_degree, mat_type;
} sparse_matrix_c;

void get_sparse_matrix_c(sparse_matrix_c *cmat, int mat_type);
void init_c(sparse_matrix_c *cmat, int nrow, int ncol, int nnz);
void build_c(sparse_matrix_c *cmat, int *rows, int *cols, int nnz);
void get_value_c(sparse_matrix_c *cmat, int i, int j, double *val);
void get_neighbors_c(sparse_matrix_c *cmat, int i, int *nbrs);
void set_value_c(sparse_matrix_c *cmat, int i, int j, double val);
void add_value_c(sparse_matrix_c *cmat, int i, int j, double val);
void matvec_c(sparse_matrix_c *cmat, double *x, double *y, int n, int m);
}
