#include <matrix.hpp>

extern"C" {

typedef struct {
    void *p;
    int nn, solver_type;
} iterative_solver_c;

typedef struct {
    void *p;
    int nn, level, pc_type;
} preconditioner_c;

void get_iterative_solver(iterative_solver_c *csolver, int solver_type);
void get_preconditioner(preconditioner_c *cpc, int pc_type);

void solver_init(iterative_solver_c *csolver, int nn, double tolerance);
void solve(iterative_solver_c *csolver, sparse_matrix_c *cmat, double *x, double *b, preconditioner_c *cpc, int n);

void preconditioner_init(preconditioner_c *cpc, sparse_matrix_c *cmat, int level);
void precondition(preconditioner_c *cpc, sparse_matrix_c *cmat, double *x, double *b, int n);

}
