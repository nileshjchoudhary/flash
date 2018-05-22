#include "fc.h"
int flash_calculation_solve_dense_linear_system_LU(int n, double *a0, int pvt[])
{
    int i, j, k;
    double d, r;

#define a(i,j)	(a0[(i) * n + (j)])

    for (i = 0; i < n; i++) {
        d = fabs(a(i,i));
        k = i;
        for (j = i + 1; j < n; j++) {
            if ((r = fabs(a(j,i))) > d) {
                d = r;
                k = j;
            }
        }

        if (d == 0.0) 
            return 0;

        pvt[i] = k;
        if (k != i) { /* exchange row i and row k */
            for (j = i; j < n; j++) {
                d = a(i,j);
                a(i,j) = a(k,j);
                a(k,j) = d;
            }
        }

        if ((d = a(i,i)) != 1.0) {
            a(i,i) = d = 1.0 / d;
            for (j = i + 1; j < n; j++)
                a(i,j) *= d;
        }

        for (j = i + 1; j < n; j++) {
            if ((d = a(j,i)) == 0.0) continue;

            for (k = i + 1; k < n; k++) a(j,k) -= d * a(i,k);
        }
    }

#undef a

    return 1;
}

/* solves L*U*X = B. */
void flash_calculation_solve_dense_linear_system_SV(int n, double *a0, 
        int pvt[], int m, double *b0)
{
    int i, j, k;
    double d;

#define a(i,j)	(a0[(i) * n + (j)])
#define b(i,j)	(b0[(i) * m + (j)])

    for (i = 0; i < n; i++) {
        k = pvt[i];
        if (k != i) {
            /* exchange row i with row k */
            for (j = 0; j < m; j++) {
                d = b(i,j);
                b(i,j) = b(k,j);
                b(k,j) = d;
            }
        }

        if ((d = a(i,i)) != 1.0) {
            for (j = 0; j < m; j++) b(i,j) *= d;
        }

        for (j = i + 1; j < n; j++) {
            if ((d = a(j,i)) == 0.0) continue;
            for (k = 0; k < m; k++) b(j,k) -= d * b(i,k);
        }
    }

    for (i = n - 2; i >= 0; i--)
        for (j = i + 1; j < n; j++) {
            d = a(i,j);
            for (k = 0; k < m; k++) b(i,k) -= d * b(j,k);
        }
#undef a
#undef b
    return;
}

int flash_calculation_solve_dense_linear_system(double *M, double *b, 
        double *x, int n)
{
    int *pvt, i, LU;

    pvt = malloc(n * sizeof(*pvt));

    LU = flash_calculation_solve_dense_linear_system_LU(n, M, pvt);

    if (LU) {
        flash_calculation_solve_dense_linear_system_SV(n, M, 
                pvt, 1, b);

        for (i = 0; i < n; i++) {
            x[i] = b[i];
        }

        free(pvt);
    }
    else {
        free(pvt);

        return LU;
    }

    return 1;
}

