#ifndef FC_SOLVER
int flash_calculation_solve_dense_linear_system_LU(int n, double *a0, int pvt[]);
void flash_calculation_solve_dense_linear_system_SV(int n, double *a0, 
        int pvt[], int m, double *b0);
int flash_calculation_solve_dense_linear_system(double *M, double *b, 
        double *x, int n);

#define FC_SOLVER
#endif
