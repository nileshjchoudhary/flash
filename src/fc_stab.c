#include "fc.h"

static int stab_itr = 0;
static double stab_solve_time = 0.;

/* ## 2. Stability Test
# The following code is designed for stability test. */

/* ### Use Wilson's equation as initial guess for K-values
# $$
# K_i = \frac{P_{c,i}}{P} \exp{5.373 (1 + \omega_i)(1 - \frac{T_{c,i}}{T}})
# $$
*/

double * flash_calculation_estimate_K(EOS *eos, double *K)
{
    int i, ncomp = eos->ncomp;
    double P, T;
    COMP *comp;

    if (K == NULL) {
        K = malloc(ncomp * sizeof(*K));
    }

    P = eos->pres;
    T = eos->temp;
    comp = eos->comp_list->comp;

    for (i = 0; i < ncomp; i++) {
        K[i] = comp[i].PC / P 
            * exp(5.37 * (1.0 + comp[i].AC) * (1.0 - comp[i].TC / T));
    }

    return K;
}

/* ### Composition Initial Guess List
# The initial guess is from the paper "General Strategy for Stability Testing and Phase-Split Calculation in Two and Three Phases" by Zhidong Li and Abbas Firoozabadi in SPE Journal, 2012.
# 
# In total, there are $N_c + 4$ sets of initial guess:
# 1. Wilson's equation: $X_i = K_i x_i$
# 2. Inverse Wilson's equation: $X_i = x_i / K_i$
# 3. $X_i = K_i^{\frac{1}{3}} x_i$
# 4. $X_i = x_i / K_i^{\frac{1}{3}}$
# 5. Pure component: 
# $$
# X_i = 
# \left\{
# \begin{array}[c]
#  0.9 \quad i = j \\
#  0.1 / (N_c - 1) \quad otherwise
# \end{array}
# \right.
# $$ for $j = 1, \cdots, N_c$
*/

double * flash_calculation_stability_analysis_initial_estimate(PHASE *phase)
{
    int ncomp = phase->ncomp, n_guess, i, j;
    double *K, Xi, *est;
    EOS *eos = phase->eos;

    est = malloc((ncomp + 4) * ncomp * sizeof(*est));
    n_guess = 0;

    K = malloc(ncomp * sizeof(*K));
    flash_calculation_estimate_K(eos, K);

    /* Wilson correlation */
    for (i = 0; i < ncomp; i++) {
        Xi = K[i] * phase->mf[i];
        *(est + n_guess * ncomp + i) = Xi;
    }
    n_guess += 1;

    /* Inverse Wilson correlation */
    for (i = 0; i < ncomp; i++) {
        Xi = phase->mf[i] / K[i];
        *(est + n_guess * ncomp + i) = Xi;
    }
    n_guess += 1;

    /* Wilson correlation power(1./3.) */
    for (i = 0; i < ncomp; i++) {
        Xi = pow(K[i], 1.0 / 3.0) * phase->mf[i];
        *(est + n_guess * ncomp + i) = Xi;
    }
    n_guess += 1;

    /* Inverse Wilson correlation power(1./3.) */
    for (i = 0; i < ncomp; i++) {
        Xi = phase->mf[i] / pow(K[i], 1.0 / 3.0); 
        *(est + n_guess * ncomp + i) = Xi;
    }
    n_guess += 1;

    /* A pure phase */
    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                Xi = 0.9;
            }
            else {
                Xi = 0.1 / (ncomp - 1);
            }
            *(est + n_guess * ncomp + j) = Xi;
        }
        n_guess += 1;
    }


    /* A hypothetical idea gas: TODO */
    free(K);

    return est;
}

/* ### Calculate compositions using X
# $$
# x_i = \frac{X_i}{\sum_j X_j} 
# $$
*/

void flash_calculation_calculate_trial_phase_composition(double *X_t, double *x, int ncomp)
{
    double Xs = 0.0;
    int i;

    for (i = 0; i < ncomp; i++) {
        Xs += X_t[i];
    }

    for (i = 0; i < ncomp; i++) {
        x[i] = X_t[i] / Xs;
    }
}

/* ### Update X using SS method
# $$
# X_{t,i} = \exp(\log{x_i} + \log{\phi_i} - \log{\phi_{t,i}})
# $$
*/

void flash_calculation_SS_method_update_X(PHASE *phase, PHASE *phase_t, double *X_t)
{
    int ncomp = phase->ncomp, i;

    for (i = 0; i < ncomp; i++) {
        X_t[i] = exp(log(phase->mf[i]) + log(phase->phi[i]) - log(phase_t->phi[i]));
    }
}

/* ### Calculate derivatives of x with respecte to X
# From $x_i = \frac{X_i}{\sum{X_j}}$, we have
# $$
# \frac{\partial x_i}{\partial X_j} = \frac{\sigma_{i,j}}{\sum{X_j}}
#     - \frac{X_i}{(\sum{X_j})^2}
# $$
# where
# $$
# \sigma_{i,j} = 
# \left\{
# \begin{array}{cc}
# 1 & i = j\\
# 0 & \text{otherwise}
# \end{array}
# \right.
# $$
*/

void flash_calculation_calculate_trial_phase_composition_derivative(double *X_t, double *dx_t, int ncomp)
{
    double sum_X = 0.0, sigma = 0.0;
    int i, j;

    for (i = 0; i < ncomp; i++) {
        sum_X += X_t[i];
    }

    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                sigma = 1.0;
            }
            else {
                sigma = 0.0;
            }

            dx_t[i * ncomp + j] = sigma / sum_X - X_t[i] / (sum_X * sum_X);
        }
    }
}

/* ### Calculate the values of equilibrium equations
# $$
#     D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)}
# $$
*/

void flash_calculation_calculate_stability_equilibrium_equation(PHASE *phase, PHASE *phase_t, 
        double *X_t, double *D)
{
    int ncomp = phase->ncomp, i;

    for (i = 0; i < ncomp; i++) {
        D[i] = log(X_t[i]) - log(phase->mf[i]) + log(phase_t->phi[i]) - log(phase->phi[i]);
    }
}

/* # ### Calculate the derivatives of equilibrium equations with respect to X
# From the equations
# $$
#     D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)},
# $$
# we have the following derivatives
# $$
# \frac{\partial D_i}{\partial X_j} = \frac{\sigma_{i,j}}{X_i} 
#                             + \frac{1}{\phi_i(x)} 
#                             \sum_k{\frac{\partial \phi_i(x)}{\partial x_k} \frac{\partial x_k}{\partial X_j}}
# $$
*/

void flash_calculation_calculate_stability_equilibrium_equation_derivative(PHASE *phase_t, double *dx_t, 
        double *X_t, double *dD)
{
    int ncomp = phase_t->ncomp, i, j, k;
    double sigma = 0.0, temp;

    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                sigma = 1.0;
            }
            else {
                sigma = 0.0;
            }

            dD[i * ncomp + j] = 1.0 / X_t[i] * sigma;

            temp = 0.0;
            for (k = 0; k < ncomp; k++) {
                temp += phase_t->dphi_dx[i * ncomp + k] * dx_t[k * ncomp + j];
            }

            dD[i * ncomp + j] += 1.0 / phase_t->phi[i] * temp;
        }
    }
}

double flash_calculation_calculate_stability_residual(PHASE *phase, PHASE *phase_t, double *X_t, double *res)
{
    int ncomp = phase->ncomp, i;
    double tol = 0.0, tmp1, tmp2;

    for (i = 0; i < ncomp; i++) {
        res[i] = log(X_t[i]) + log(phase_t->phi[i]) - log(phase->mf[i]) - log(phase->phi[i]);

        tmp1 = res[i] * res[i];
        tmp2 = log(phase->mf[i]) + log(phase->phi[i]);
        tmp2 = tmp2 * tmp2;

        if ((tmp1 / tmp2) > tol) {
            tol = tmp1 / tmp2;
        }
    }

    return tol;
}


/* ### Update X using QNSS method
# Solve $\delta X$ through the following equation
# $$
# J \delta X = - D,
# $$
# where $D_i$ is the value of equilibrium equation, 
# $J$ is Jacobian matrix and $J_{i,j} = \frac{\partial D_i}{\partial X_j}$,
# and update X by 
# $$
# X = X + \delta X
# $$
*/

void flash_calculation_QNSS_method_update_X(double *dD, double *D, double *X_t, int ncomp)
{
    int i;
    double *x;

    x = malloc(ncomp * sizeof(*x));

    for (i = 0; i < ncomp; i++) {
        D[i] = - D[i];
    }

    flash_calculation_solve_dense_linear_system(dD, D, x, ncomp);

    for (i = 0; i < ncomp; i++) {
        X_t[i] += x[i];
    }

    free(x);
}

/* ### Check Stability using X
# 1. If $\sum_i{(\log{\frac{X_i}{z_i}})^2} < \epsilon$, the solution is trivial, try next initial guess.
# 2. If $\sum_i X_i < 1.0$, the phase is stable; otherwise, it is unstable.
*/

int flash_calculation_check_stability(double *X_t, double *z, int ncomp)
{
    int i;
    double sum_Y = 0.0, sum_K = 0.0;

    for (i = 0; i < ncomp; i++) { 
        sum_Y += X_t[i];

        sum_K += log(X_t[i] / z[i]) * log(X_t[i] / z[i]);
    }

    if (sum_K < 1e-2) 
        return -1;

    if (sum_Y < 1.0 + 1e-8) {
        return 1;
    }
    else {
        return 0;
    }

    return -1;
}

/* ### QNSS method to solve stability analysis
# The tangent-plane criterion of stability analysis of a phase with 
#         compositions z results in solving the following set of equations 
#         (Nghiem and Li, 1984):
# $$
#             D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)}
# $$
# where,
# $$
#             x_i = X_i / \sum{X_j}
# $$
# for the primary variables X. 
#         The phase will be stable if (1) $\sum{X_i} < 1$
#                                     (2) $\sum{X_i} = 1$ and $X_i \neq z_i$
#         otherwise the phase is unstable.
#         Equations D are solved using the SS method first and then QNSS 
#         method.
#         To solve the above equations, several sets of initial guesses 
#         should be used to avoid trivial solutions.
*/

int flash_calculation_stability_analysis_QNSS(PHASE *phase, double *K, double tol)
{
    int ncomp = phase->ncomp, i, j, n_guess, itr = 0;
    double *D, *dD, *res, *est;
    double *x_t, *dx_t, *X_t, error;
    PHASE *phase_t;
    int system_status = -1;
    double solve_time;

    solve_time = flash_calculation_get_time(NULL);

    D = malloc(ncomp * sizeof(*D));
    dD = malloc(ncomp * ncomp * sizeof(*dD));
    res = malloc(ncomp * sizeof(*res));

    n_guess = ncomp + 4;
    est = flash_calculation_stability_analysis_initial_estimate(phase);

    x_t = malloc(ncomp * sizeof(*x_t));
    X_t = malloc(ncomp * sizeof(*X_t));
    dx_t = malloc(ncomp * ncomp * sizeof(*dx_t));

    phase_t = flash_calculation_phase_new(phase->eos, x_t);

    flash_calculation_compute_phase_parameter(phase);
    flash_calculation_calculate_compressibility_factor(phase);
    flash_calculation_calculate_fugacity(phase);

    for (i = 0; i < n_guess; i++) {
        for (j = 0; j < ncomp; j++) {
            X_t[j] = est[i * ncomp + j];
        }

        itr = 0;

        while(1) {
            /* Calculate trial phase composition */
            flash_calculation_calculate_trial_phase_composition(X_t, phase_t->mf, ncomp);

            /* Calculation trial phase compostion derivative */
            flash_calculation_calculate_trial_phase_composition_derivative(X_t, dx_t, ncomp);

            /* Calculate the compressibility factor and 
               fugacity coefficient for the trial phase 
               with compositions x_t */
            flash_calculation_compute_phase_parameter(phase_t);
            flash_calculation_calculate_compressibility_factor(phase_t);
            flash_calculation_calculate_fugacity(phase_t);

            /* Calculate the residual */
            error = flash_calculation_calculate_stability_residual(phase, phase_t, X_t, res);

            /* Check if stop */
            if (error < tol)
                break;

            /* Update X_t */
            if (error > 1e-5) {
                flash_calculation_SS_method_update_X(phase, phase_t, X_t);
            }
            else {
                /* Calculate the equilibrim equation and its derivative */
                flash_calculation_calculate_stability_equilibrium_equation(phase, phase_t, X_t, D);
                flash_calculation_calculate_stability_equilibrium_equation_derivative(phase_t, dx_t, X_t, dD);

                /* Update X_t by X_t += - dD^(-1) * D */
                flash_calculation_QNSS_method_update_X(dD, D, X_t, ncomp);
            }

            /* Maximum itrations */
            itr += 1;
            if (itr > 1000) {
                //printf("##### WARNING: Stability_analysis_QNSS reach the maximum iterations!\n");
                break;
            }
        }

        /* Check stability based on Sum(X_t) */
        system_status = flash_calculation_check_stability(X_t, phase->mf, ncomp);

        /* If the solution is trivial, try the next initial guess; 
           otherwise break 
           if system_status is 'Unstable' or system_status is 'Stable':
           break */
        if (system_status != -1) {
            break;
        }
    }

    stab_itr += itr;

    /* if K is not None, we output K = X_t / z */
    if (K != NULL) {
        for (i = 0; i < ncomp; i++) {
            K[i] = X_t[i] / phase->mf[i];
        }
    }

    if (system_status == -1) {
        system_status = 1;
    }


    free(x_t);
    free(X_t);
    free(dx_t);
    free(D);
    free(dD);
    free(res);
    free(est);

    flash_calculation_phase_free(&phase_t);

    solve_time = flash_calculation_get_time(NULL) - solve_time;
    stab_solve_time += solve_time;

    return system_status;
}


/* ## Draw Stability Test Map
# For a feed composition, this function will draw a stability test map which 
    shows the results from stability tests at the given temperature range and 
    pressure range. For the example, the following figure shows the results of 
    stability tests for the temperature in [200 K, 340 K] and the pressure in 
    [20 atm, 80 atm].
*/

void flash_calculation_output_stability_analysis_map(STABILITY_MAP *sm, 
        double *comp_X, int ncomp, char *output_name)
{
    int i, j;
    char file_unstable[100], file_liquid[100], file_vapor[100], file_combination[100];
    FILE *fp;

    sprintf(file_unstable, "%s-unstable.csv", output_name);
    fp = fopen(file_unstable, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure\n");

    for (i = 0; i < sm->n_unstable; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf\n", sm->unstable_temp[i], sm->unstable_pres[i]);
    }
    fclose(fp);

    sprintf(file_liquid, "%s-liquid.csv", output_name);
    fp = fopen(file_liquid, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure\n");

    for (i = 0; i < sm->n_liquid; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf\n", sm->liquid_temp[i], sm->liquid_pres[i]);
    }
    fclose(fp);

    sprintf(file_vapor, "%s-vapor.csv", output_name);
    fp = fopen(file_vapor, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure\n");

    for (i = 0; i < sm->n_vapor; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf\n", sm->vapor_temp[i], sm->vapor_pres[i]);
    }
    fclose(fp);

    sprintf(file_combination, "%s-stability-combination.csv", output_name);
    fp = fopen(file_combination, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure,Unstable,Liquid,Vapor\n");
    for (i = 0; i < sm->n_unstable; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf,1,0,0\n", sm->unstable_temp[i], sm->unstable_pres[i]);
    }

    for (i = 0; i < sm->n_liquid; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf,0,1,0\n", sm->liquid_temp[i], sm->liquid_pres[i]);
    }

    for (i = 0; i < sm->n_vapor; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf,0,0,1\n", sm->vapor_temp[i], sm->vapor_pres[i]);
    }
    fclose(fp);
}

void flash_calculation_output_stability_analysis_map_PM(STABILITY_PM_MAP *sm, 
        double *comp_X, int ncomp, char *output_name)
{
    int i, j;
    char file_unstable[100], file_liquid[100], file_vapor[100], 
         file_combination[100];
    FILE *fp;

    sprintf(file_unstable, "%s-unstable-PM.csv", output_name);
    fp = fopen(file_unstable, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure\n");

    for (i = 0; i < sm->n_unstable; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->unstable_x[i][j]);
            }
        }
        fprintf(fp, "%lf\n", sm->unstable_pres[i]);
    }
    fclose(fp);

    sprintf(file_liquid, "%s-liquid-PM.csv", output_name);
    fp = fopen(file_liquid, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure\n");

    for (i = 0; i < sm->n_liquid; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->liquid_x[i][j]);
            }
        }
        fprintf(fp, "%lf\n", sm->liquid_pres[i]);
    }
    fclose(fp);

    sprintf(file_vapor, "%s-vapor-PM.csv", output_name);
    fp = fopen(file_vapor, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure\n");

    for (i = 0; i < sm->n_vapor; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->vapor_x[i][j]);
            }
        }
        fprintf(fp, "%lf\n", sm->vapor_pres[i]);
    }
    fclose(fp);

    sprintf(file_combination, "%s-stability-combination-PM.csv", output_name);
    fp = fopen(file_combination, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure,Unstable,Liquid,Vapor\n");
    for (i = 0; i < sm->n_unstable; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->unstable_x[i][j]);
            }
        }
        fprintf(fp, "%lf,1,0,0\n", sm->unstable_pres[i]);
    }

    for (i = 0; i < sm->n_liquid; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->liquid_x[i][j]);
            }
        }
        fprintf(fp, "%lf,0,1,0\n", sm->liquid_pres[i]);
    }

    for (i = 0; i < sm->n_vapor; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->vapor_x[i][j]);
            }
        }
        fprintf(fp, "%lf,0,0,1\n", sm->vapor_pres[i]);
    }
    fclose(fp);
}

STABILITY_MAP * flash_calculation_draw_stability_analysis_map(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double dT, double dP, FLASH_STAB_ANN *fsa, char *output_name)
{
    int n_pres, n_temp, i, j, ncomp = comp_list->ncomp;
    double *pres_list, *temp_list;
    EOS *eos;
    PHASE *phase;
    STABILITY_MAP *sm;
    int status;

    n_pres = (int)((P_max - P_min) / dP);
    n_temp = (int)((T_max - T_min) / dT);

    pres_list = malloc(n_pres * sizeof(*pres_list));
    temp_list = malloc(n_temp * sizeof(*temp_list));

    for (i = 0; i < n_pres; i++) {
        pres_list[i] = P_min + i * dP;
    }

    for (i = 0; i < n_temp; i++) {
        temp_list[i] = T_min + i * dT;
    }

    sm = malloc(sizeof(*sm));
    sm->n_unstable = 0;
    sm->unstable_pres = malloc(n_pres * n_temp * sizeof(*(sm->unstable_pres)));
    sm->unstable_temp = malloc(n_pres * n_temp * sizeof(*(sm->unstable_temp)));
    sm->n_liquid = 0;
    sm->liquid_pres = malloc(n_pres * n_temp * sizeof(*(sm->liquid_pres)));
    sm->liquid_temp = malloc(n_pres * n_temp * sizeof(*(sm->liquid_temp)));
    sm->n_vapor = 0;
    sm->vapor_pres = malloc(n_pres * n_temp * sizeof(*(sm->vapor_pres)));
    sm->vapor_temp = malloc(n_pres * n_temp * sizeof(*(sm->vapor_temp)));

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);
    phase = flash_calculation_phase_new(eos, comp_X);

    for (i = 0; i < n_pres; i++) {
        for (j = 0; j < n_temp; j++) {
            int flag = 0;

            eos->pres = pres_list[i];
            eos->temp = temp_list[j];

            if (fsa != NULL) {
                double *input;
                int n, k;

                n = ncomp + 2;
                input = malloc(n * sizeof(*input));

                for (k = 0; k < ncomp; k++) {
                    input[k] = comp_X[k];
                }
                input[ncomp] = temp_list[j];
                input[ncomp + 1] = pres_list[j];

                flag = flash_calculation_stab_ann_predict(fsa, input, n, &status);
            }

            if (!flag) {
                status = flash_calculation_stability_analysis_QNSS(phase, NULL, 1e-10);
            }

            if (status == 1) {
                if (phase->phase_no == 0) {
                    int k;

                    k = sm->n_liquid;
                    sm->liquid_temp[k] = temp_list[j];
                    sm->liquid_pres[k] = pres_list[i];

                    sm->n_liquid += 1;
                }
                else {
                    int k;

                    k = sm->n_vapor;
                    sm->vapor_temp[k] = temp_list[j];
                    sm->vapor_pres[k] = pres_list[i];

                    sm->n_vapor += 1;
                }
            }
            else if (status == 0) {
                int k;

                k = sm->n_unstable;
                sm->unstable_temp[k] = temp_list[j];
                sm->unstable_pres[k] = pres_list[i];

                sm->n_unstable += 1;
            }
        }
    }

    if (output_name != NULL) {
        flash_calculation_output_stability_analysis_map(sm, 
                comp_X, ncomp, output_name);
    }

    free(pres_list);
    free(temp_list);
    flash_calculation_phase_free(&phase);

    free(eos);

    return sm;
}

STABILITY_PM_MAP * flash_calculation_draw_stability_analysis_map_PM(COMP_LIST *comp_list, 
        double *comp_X, double T, double P_min, double P_max, double dP, int selected_component, 
        double dxx, FLASH_STAB_ANN *fsa, char *output_name)
{
    int n_pres, n_x, i, j, k, ncomp = comp_list->ncomp;
    double *pres_list, *x_list, *x, sum_no_selected;
    EOS *eos;
    PHASE *phase;
    STABILITY_PM_MAP *sm;
    int status;

    x = malloc(ncomp * sizeof(*x));
    n_pres = (int)((P_max - P_min) / dP);
    pres_list = malloc(n_pres * sizeof(*pres_list));
    for (i = 0; i < n_pres; i++) {
        pres_list[i] = P_min + i * dP;
    }

    n_x = (int)((1.0 - 0.001) / dxx) + 1;
    x_list = malloc(n_x * sizeof(*x_list));
    for (i = 0; i < n_x - 1; i++) {
        x_list[i] = 0.001 + i * dxx;
    }
    x_list[n_x - 1] = 0.999;

    sm = malloc(sizeof(*sm));
    sm->n_unstable = 0;
    sm->unstable_pres = malloc(n_pres * n_x * sizeof(*(sm->unstable_pres)));
    sm->unstable_x = malloc(n_pres * n_x * sizeof(*(sm->unstable_x)));
    sm->n_liquid = 0;
    sm->liquid_pres = malloc(n_pres * n_x * sizeof(*(sm->liquid_pres)));
    sm->liquid_x = malloc(n_pres * n_x * sizeof(*(sm->liquid_x)));
    sm->n_vapor = 0;
    sm->vapor_pres = malloc(n_pres * n_x * sizeof(*(sm->vapor_pres)));
    sm->vapor_x = malloc(n_pres * n_x * sizeof(*(sm->vapor_x)));

    for (i = 0; i < n_pres * n_x; i++) {
        sm->unstable_x[i] = malloc(ncomp * sizeof(double));
        sm->liquid_x[i] = malloc(ncomp * sizeof(double));
        sm->vapor_x[i] = malloc(ncomp * sizeof(double));
    }

    eos = flash_calculation_EOS_new(comp_list, 0.0, T, 0);
    phase = flash_calculation_phase_new(eos, x);

    sum_no_selected = 0.0;
    for (i = 0; i < ncomp; i++) {
        if (i != selected_component) {
            sum_no_selected += comp_X[i];
        }
    }

    for (i = 0; i < n_pres; i++) {
        for (j = 0; j < n_x; j++) {
            int flag = 0;

            for (k = 0; k < ncomp; k++) {
                if (k == selected_component) {
                    x[k] = x_list[j];
                }
                else {
                    x[k] = (1.0 - x_list[j]) * comp_X[k] / sum_no_selected;
                }
            }

            eos->pres = pres_list[i];
            eos->temp = T;

            if (fsa != NULL) {
                double *input;
                int n, k;

                n = ncomp + 1;
                input = malloc(n * sizeof(*input));

                for (k = 0; k < ncomp; k++) {
                    input[k] = x[k];
                }
                input[ncomp] = pres_list[i];

                flag = flash_calculation_stab_ann_predict(fsa, input, n, &status);
            }

            if (!flag) {
                status = flash_calculation_stability_analysis_QNSS(phase, NULL, 1e-10);
            }


            if (status == 1) {
                if (phase->phase_no == 0) {
                    int l;

                    l = sm->n_liquid;
                    sm->liquid_pres[l] = pres_list[i]; 

                    for (k = 0; k < ncomp; k++) {
                        sm->liquid_x[l][k] = x[k];
                    }

                    sm->n_liquid += 1;
                }
                else {
                    int l;

                    l = sm->n_vapor;
                    sm->vapor_pres[l] = pres_list[i];

                    for (k = 0; k < ncomp; k++) {
                        sm->vapor_x[l][k] = x[k];
                    }

                    sm->n_vapor += 1;
                }
            }
            else if (status == 0) {
                int l;

                l = sm->n_unstable;
                sm->unstable_pres[l] = pres_list[i];

                for (k = 0; k < ncomp; k++) {
                    sm->unstable_x[l][k] = x[k];
                }

                sm->n_unstable += 1;
            }
        }
    }

    if (output_name != NULL) {
        flash_calculation_output_stability_analysis_map_PM(sm, 
                NULL, ncomp, output_name);
    }

    free(pres_list);
    free(x_list);
    free(x);
    flash_calculation_phase_free(&phase);

    free(eos);

    return sm;
}

void flash_calculation_stability_map_free(STABILITY_MAP **sm)
{
    STABILITY_MAP *sm0 = *sm;

    free(sm0->unstable_pres);
    free(sm0->unstable_temp);

    free(sm0->liquid_pres);
    free(sm0->liquid_temp);

    free(sm0->vapor_pres);
    free(sm0->vapor_temp);

    free(*sm);
}

void flash_calculation_stability_PM_map_free(STABILITY_PM_MAP **sm)
{
    int i;
    STABILITY_PM_MAP *sm0 = *sm;

    for (i = 0; i < sm0->n_unstable; i++) {
        free(sm0->unstable_x[i]);
    }
    free(sm0->unstable_pres);
    free(sm0->unstable_x);

    for (i = 0; i < sm0->n_unstable; i++) {
        free(sm0->liquid_x[i]);
    }
    free(sm0->liquid_pres);
    free(sm0->liquid_x);

    for (i = 0; i < sm0->n_unstable; i++) {
        free(sm0->vapor_x[i]);
    }
    free(sm0->vapor_pres);
    free(sm0->vapor_x);

    free(*sm);
}

double flash_calculation_stability_time_cost()
{
    double stab_solve_time0 = stab_solve_time;

    MPI_Allreduce(&stab_solve_time0, &stab_solve_time, 
            1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return stab_solve_time;
}

int flash_calculation_stability_iteration_number()
{
    int stab_itr0 = stab_itr;

    MPI_Allreduce(&stab_itr0, &stab_itr, 1, MPI_INT, 
            MPI_SUM, MPI_COMM_WORLD);

    return stab_itr;
}
