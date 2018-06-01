#include "fc.h"

static int verb = 0;
static int split_failure = 0;
static int split_itr = 0;
static double split_solve_time = 0.;
static double split_pred_time = 0.;

/* ## 3. Two-phase Flash Calculation
# The following code is designed for two-phase flash calculations.
*/

/* ### Calculate compositions for liquid and vapour phases
# $$
# x_{l,i} = \frac{z_i}{1 + (K_i - 1) F_v} \\
# x_{v,i} = \frac{K_i z_i}{1 + (K_i - 1) F_v}
# $$
# where $z_i$ is the feed composition, $F_v$ is mole fraction of vapour phase.
*/

void flash_calculation_calculate_composition(double *K, double *z, double F_v, 
        double *x_l, double *x_v, int ncomp)
{
    int i;

    for (i = 0; i < ncomp; i++) {
        x_l[i] = z[i] / (1.0 + (K[i] - 1.0) * F_v);
        x_v[i] = K[i] * z[i] / (1.0 + (K[i] - 1.0) * F_v);
    }
}

/* ### Calculate derivatives of compositions for liquid and vapour phases with respect to $K_i$
# $$
# \frac{\partial x_{l,i}}{\partial K_i} = - \frac{z_i F_v}{(1 + (K_i - 1) F_v)^2}  \\
# \frac{\partial x_{v,i}}{\partial K_i} = \frac{z_i}{1 + (K_i - 1) F_v} - \frac{z_i K_i F_v}{(1 + (K_i - 1) F_v)^2}
# $$
*/

void flash_calculation_calculate_composition_derivative(double *K, double *z, double F_v, 
        double *dx_l, double *dx_v, int ncomp)
{
    int i;
    double temp;

    for (i = 0; i < ncomp; i++) {
        temp = 1.0 + (K[i] - 1.0) * F_v;
        dx_l[i] = - z[i] / (temp * temp) * F_v;
        dx_v[i] = z[i] / temp - z[i] * K[i] / (temp * temp) * F_v;
    }
}

/* ### Calculate equilibrium equation
# $$
# G_i = \log{x_{v,i}} + \log{\phi_{v,i}(x_v)} 
#             - \log{x_{l,i}} - \log{\phi_{l,i}(x_l)}
# $$
*/

double flash_calculation_calculate_equilibrium_equation(PHASE *phase_L, PHASE *phase_V, double *G)
{
    int ncomp = phase_L->ncomp, i;
    double error = 0.0;

    for (i = 0; i < ncomp; i++) {
        G[i] = log(phase_V->mf[i]) + log(phase_V->phi[i]);
        G[i] += - log(phase_L->mf[i]) - log(phase_L->phi[i]);

        error += G[i] * G[i];
    }

    return error;
}

/* ### Calculate derivatives of equilibrium equations with respect to K
# $$
# \frac{\partial G_i}{\partial K_j} = \frac{1}{x_{v,i}} \frac{\partial x_{v,i}}{\partial K_i} \sigma_{i,j} 
#                 + \frac{1}{\phi_{v,i}} \frac{\partial \phi_{v,i}}{\partial x_{v,j}} \frac{\partial x_{v,j}}{\partial K_j}
#                - \frac{1}{x_{l,i}} \frac{\partial x_{l,i}}{\partial K_i} \sigma_{i,j} 
#                - \frac{1}{\phi_{l,i}} \frac{\partial \phi_{l,i}}{\partial x_{l,j}} \frac{\partial x_{l,j}}{\partial K_j}     
# $$
*/

void flash_calculation_calculate_equilibrium_equation_derivative(PHASE *phase_L, PHASE *phase_V, 
        double *dx_l, double *dx_v, double *dG)
{
    int ncomp = phase_L->ncomp, i, j;

    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                dG[i * ncomp + j] = 1.0 / phase_V->mf[i] * dx_v[i] - 1.0 / phase_L->mf[i] * dx_l[i];
            }
            else {
                dG[i * ncomp + j] = 0.0;
            }

            dG[i * ncomp + j] += 1.0 / phase_V->phi[i] * phase_V->dphi_dx[i * ncomp + j] * dx_v[j]
                - 1.0 / phase_L->phi[i] * phase_L->dphi_dx[i * ncomp + j] * dx_l[j];
        }
    }
}

/* ### Update K using QNSS method
# Solve
# $$
# dG \delta K = - G,
# $$
# and update $K$ by
# $$
# K = K + \delta K
# $$
*/

void flash_calculation_QNSS_method_update_K(double *dG, double *G, double *K, int ncomp)
{
    int i;
    double *x;

    x = malloc(ncomp * sizeof(*x));

    for (i = 0; i < ncomp; i++) {
        G[i] = -G[i];
    }

    flash_calculation_solve_dense_linear_system(dG, G, x, ncomp);

    for (i = 0; i < ncomp; i++) {
        K[i] += x[i];
    }

    free(x);
}

/* ### Calculate the value of the Rachford-Rice equation
# $$
# V = \sum_i{\frac{z_i (K_i - 1)}{1 + (K_i - 1) F_v}}
# $$
*/

double flash_calculation_calculate_RachfordRice_equation_value(double *K, double *z, double n_V, int ncomp)
{
    int i;
    double value = 0.0;

    for (i = 0; i < ncomp; i++) {
        value += z[i] * (K[i] - 1.0) / (1.0 + (K[i] - 1.0) * n_V);
    }

    return value;
}

/* ### Calculate the derivative of Rachford-Rice equation
# $$
# \frac{\partial V}{\partial F_v} = - \sum_i{\frac{z_i (K_i - 1)^2}{(1 + (K_i - 1)F_v)^2}}
# $$
*/

double flash_calculation_calculate_RachfordRice_equation_derivative(double *K, double *z, double n_V, int ncomp)
{
    int i;
    double value = 0.0;

    for (i = 0; i < ncomp; i++) {
        value += - z[i] * (K[i] - 1.0) * (K[i] - 1.0) 
            / (1.0 + (K[i] - 1.0) * n_V) 
            / (1.0 + (K[i] - 1.0) * n_V);
    }

    return value;
}

/* ### Solve the Rachford-Rice equation
# Solve $F_v$ using an iterative method:
# 1. Solve $\delta F_v$
# $$
# \frac{\partial V}{\partial F_v} \delta F_v = - V 
# $$
# 2. Update $F_v$
# $$
# F_v = F_v + \delta F_v
# $$
# 3. Re-calculate $V$ and $\frac{\partial V}{\partial F_v}$
*/

double flash_calculation_solve_RachfordRice_equation(double *K, double *z, double n_V0, int ncomp)
{
    int itr = 0;
    double n_V = n_V0, F, J, d;

    while(1) {
        F = flash_calculation_calculate_RachfordRice_equation_value(K, z, n_V, ncomp);

        if (fabs(F) < 1e-10)
            break;

        J = flash_calculation_calculate_RachfordRice_equation_derivative(K, z, n_V, ncomp);

        d = - F / J;
        n_V += d;

        if (n_V < 0.0) {       
            n_V -= d;
            n_V *= 0.5;
		}

        if (n_V > 1.0) {
            n_V -= d;
            n_V = (1.0 + n_V) * 0.5;
		}
            
        itr += 1;
        if (itr > 100)
            break;
	}       
            
    return n_V;
}

/* ### Update K using SS method
# $$
# K_i = K_i \frac{\phi_{l,i}}{\phi_{v,i}}
# $$
*/
void flash_calculation_SS_method_update_K(double *fug_L, double *fug_V, double *K, int ncomp)
{
	int i;

    for (i = 0; i < ncomp; i++) {
        K[i] = K[i] * fug_L[i] / fug_V[i];
	}
}


/*
    At range [0, 1.0] find a value which makes Rachford-Rich is zero
    using bisection method
*/
double flash_calculation_two_phase_flash_calculation_calculate_initial_F(double *K, double *z, int ncomp) 
{
	double tol = 1e-3, F_min, F_max, 
           V_min, V_max, F_mid, V_mid;
	int itr;
   
    F_min = 0.0;
    F_max = 1.0;
    V_min = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_min, ncomp);
    V_max = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_max, ncomp);
    F_mid = 0.0;
    itr = 0;
    
    while(1) {
        /* calculation the value at bisection position */ 
        F_mid = (F_min + F_max) * 0.5;
        V_mid = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_mid, ncomp);
        
        /* if Value is zero, break */
        if (fabs(V_mid) < 1e-10)
            break;
            
        /* # If V_mid has different signal from V_min,
           #  set F_mid as F_max and V_mid as V_max */
        if (V_mid * V_min < 0.0) {
            V_max = V_mid;
            F_max = F_mid;
		}
        
        /* # If V_mid has different signal from V_max,
           #  set F_mid as F_min and V_mid as V_min */
        if (V_mid * V_max < 0.0) {
            V_min = V_mid;
            F_min = F_mid;
		}
            
        if (fabs(F_min - F_max) < tol) 
            break;
        
        itr += 1;
        if (itr > 100) 
            break;
	}
    
    return F_mid;
}

/* ### QNSS method for Two-phase flash calculation
# Two-phase flahs calculation requires the solution of the following 
# equilibrium and material balance equaions:
# Equilibrium equation:
# $$
#             G_i = \log{x_{v,i}} + \log{\phi_{v,i}(x_v)} 
#             - \log{x_{l,i}} - \log{\phi_{l,i}(x_l)}
# $$
# Material balance equation:
# $$
# V = \sum_i{\frac{z_i (K_i - 1)}{1 + (K_i - 1) F_v}}
# $$
# We take $K_i$ and $F_v$ as the primary variables.
# 
# K-values are solved vis the equilibrium equations, then solve the material balance equation to get the $F_v$. Repeat the process until converge. 
*/
/*
        Two-phase flahs calculation requires the solution of the following 
        equilibrium and material balance equaions:
        Equilibrium equation:
            G[i] = ln(x_v[i]) + ln(phi_v[i]) - ln(x_l[i]) - ln(phi_l[i]) 
            for i = 1,...,Nc
        Material balance equation:
            G[Nc + 1] = Sum_k((z[k] * (K[k] - 1.0)) / (1.0 + F_v * (K[k] - 1.0)))
        We take K[i] and F_v as the primary variables.
        
*/
double flash_calculation_two_phase_flash_Calculation_QNSS(EOS *eos, double *z, 
        double *K, double Fv, double tol)
{
	int i, ncomp = eos->ncomp, itr;
	double F_v, sum_K, error, *x_l, *x_v;
	double *G, *dG, *dx_l, *dx_v, *K0;
    PHASE *phase_L, *phase_V;
    
    /* Initial estimate K */
	K0 = malloc(ncomp * sizeof(*K0));
    if (Fv <= 0.0) {
        flash_calculation_estimate_K(eos, K0);
        F_v = 0.5;

#if 0
        printf("KI: \n");
        for (i = 0; i < ncomp; i++) {
            printf("%e ", K0[i]);
        }
        printf("\n");
#endif
	}
    else {
        sum_K = 0.0;
        for (i = 0; i < ncomp; i++) {
            sum_K += log(K[i]) * log(K[i]);
		}

        if (sum_K < 1e-5) {
            flash_calculation_estimate_K(eos, K0);
            F_v = 0.5;
		}
        else {
			for (i = 0; i < ncomp; i++) 
				K0[i] = K[i];
            F_v = Fv;
		}
	}

#if 0
    {
        double *KK;

        KK = malloc(ncomp * sizeof(*KK));
        flash_calculation_estimate_K(eos, KK);

        printf("Estimate K:\n");
        for (i = 0; i < ncomp; i++) {
            printf("%e ", KK[i]);
        }
        printf("\n");

        free(KK);
    }
#endif

    /* Initial compositions x_v and x_l */
	x_l = malloc(ncomp * sizeof(*x_l));
	x_v = malloc(ncomp * sizeof(*x_v));
    flash_calculation_calculate_composition(K0, z, F_v, x_l, x_v, ncomp);
   
    phase_L = flash_calculation_phase_new(eos, x_l);
    phase_V = flash_calculation_phase_new(eos, x_v);
    
    G = malloc(ncomp * sizeof(*G));
    dG = malloc(ncomp * ncomp * sizeof(*dG));
    dx_l = malloc(ncomp * sizeof(*dx_l));
    dx_v = malloc(ncomp * sizeof(*dx_v));
    
    itr = 0;

    while(1) {
        /* Calculate liquid phase fugacity */
        flash_calculation_compute_phase_parameter(phase_L);
        flash_calculation_calculate_compressibility_factor(phase_L);
        flash_calculation_calculate_fugacity(phase_L);
        
        /* Calculate vapour phase fugacity */
        flash_calculation_compute_phase_parameter(phase_V);
        flash_calculation_calculate_compressibility_factor(phase_V);
        flash_calculation_calculate_fugacity(phase_V);
        
        /* Calculate G error */
        error = flash_calculation_calculate_equilibrium_equation(phase_L, phase_V, G);
        
        /* Check convergence */
        if (error < tol)
            break;
        
        if (error > 1e-5) {
            flash_calculation_SS_method_update_K(phase_L->fug, phase_V->fug, K0, ncomp);
		}
        else {
            /* Calculate the derivatives */
            flash_calculation_calculate_composition_derivative(K0, z, F_v, dx_l, dx_v, ncomp);
            flash_calculation_calculate_equilibrium_equation_derivative(phase_L, phase_V, dx_l, dx_v, dG);
        
            /* Update K */
            flash_calculation_QNSS_method_update_K(dG, G, K0, ncomp);
		}
        
        /* ## Solve Rachford-Rice equation and get F_v */
        F_v = flash_calculation_solve_RachfordRice_equation(K0, z, F_v, ncomp);
        
        /* ## Calculate compositions */
        flash_calculation_calculate_composition(K0, z, F_v, x_l, x_v, ncomp);
        
        itr += 1;
        if (itr > 1000) {
            //printf("##### WARNING: two-phase flash calculation reach maximum iterations!\n");
            break;
		}
	}

    if (verb) {
        printf("Pres: %e, Temp: %e, Itr: %d\n", eos->pres, eos->temp, itr);
    }

    if (fabs(F_v) < 1e-5 || fabs(F_v - 1.0) < 1e-5) {
        split_failure ++;
    }
    split_itr += itr;

    flash_calculation_calculate_phase_density(phase_V);
    flash_calculation_calculate_phase_density(phase_L);

    if (phase_V->density > phase_L->density) {
        F_v = 1.0 - F_v;
        for (i = 0; i < ncomp; i++) {
            K[i] = 1.0 / K0[i];
        }
    }
    else {
        for (i = 0; i < ncomp; i++) {
            K[i] = K0[i];
        }
    }
	
	free(x_l);
	free(x_v);
	free(G);
	free(dG);
	free(dx_l);
	free(dx_v);
	free(K0);

    flash_calculation_phase_free(&phase_L);
    flash_calculation_phase_free(&phase_V);
    /* printf("##### Two-phase flash calculation iterations: %d" %itr); */
    
    return F_v;
}




/* ## Draw Two-phase Flash Calculation Map
# For a feed composition, this function will draw a two-phase flash calculation 
map which shows mole fraction of vapour phase at the given temperature range 
and pressure range. For the example, the following figure shows the results 
of two-phase flash calculation for the temperature in [200 K, 340 K] and the 
pressure in [20 atm, 80 atm].
*/

void flash_calculation_output_split_calculation_map(SPLIT_MAP *sm, 
        double *comp_X, int ncomp, int filter, char *output_name)
{
    int i, j;
    char file_name[100];
    FILE *fp;

    sprintf(file_name, "%s-split-calculation.csv", output_name);
    fp = fopen(file_name, "a");
    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure,Fv");
    for (j = 0; j < ncomp; j++) {
        fprintf(fp, ",K_%d", j);
    }
    fprintf(fp, "\n");

    for (i = 0; i < sm->n; i++) {
        if (filter && (fabs(sm->F[i]) < 1e-5 || fabs(sm->F[i] - 1.0) < 1e-5)) {
            continue;
        }

        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", sm->x[i][j]);
            }
        }

        fprintf(fp, "%lf,%lf,%lf", sm->temp[i], sm->pres[i], sm->F[i]);

        for (j = 0; j < ncomp; j++) {
            fprintf(fp, ",%lf", sm->K[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void flash_calculation_output_split_calculation_map_PM(SPLIT_PM_MAP *sm, 
        double *comp_X, int ncomp, int filter, char *output_name)
{
    int i, j;
    char file_name[100];
    FILE *fp;

    sprintf(file_name, "%s-split-calculation-PM.csv", output_name);
    fp = fopen(file_name, "a");
    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure,Fv");
    for (j = 0; j < ncomp; j++) {
        fprintf(fp, ",K_%d", j);
    }
    fprintf(fp, "\n");

    for (i = 0; i < sm->n; i++) {
        if (filter && (fabs(sm->F[i]) < 1e-5 || fabs(sm->F[i] - 1.0) < 1e-5)) {
            continue;
        }

        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", sm->x[i][j]);
            }
        }

        fprintf(fp, "%lf,%lf", sm->pres[i], sm->F[i]);

        for (j = 0; j < ncomp; j++) {
            fprintf(fp, ",%lf", sm->K[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

SPLIT_MAP * flash_calculation_draw_split_calculation_map(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double dT, double dP, FLASH_SPLIT_ANN *fsa, char *output_name)
{
    int i, j, k, ncomp = comp_list->ncomp;
    int n_pres, n_temp;
    double *pres_list, *temp_list;
    SPLIT_MAP *sm;
    int status;
    double Fv, Fv0, *K0, *K, solve_time;
    EOS *eos;
    PHASE *phase;

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
    sm->n = n_pres * n_temp;
    sm->temp = malloc(n_pres * n_temp * sizeof(*(sm->temp)));
    sm->pres = malloc(n_pres * n_temp * sizeof(*(sm->pres)));
    sm->F = malloc(n_pres * n_temp * sizeof(*(sm->F)));
    sm->K = malloc(n_pres * n_temp * sizeof(*(sm->K)));
    sm->x = malloc(n_pres * n_temp * sizeof(*(sm->x)));
    for (i = 0; i < n_pres * n_temp; i++) {
        sm->K[i] = malloc(ncomp * sizeof(*(sm->K[i])));
        sm->x[i] = malloc(ncomp * sizeof(*(sm->x[i])));
    }

    K0 = malloc(ncomp * sizeof(*K0));
    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);
    phase = flash_calculation_phase_new(eos, comp_X);

    for (i = 0; i < n_pres; i++) {
        K = NULL;
        Fv = Fv0 = 0.0;

        for (j = 0; j < n_temp; j++) {
            eos->pres = pres_list[i];
            eos->temp = temp_list[j];

            sm->pres[i * n_temp + j] = pres_list[i];
            sm->temp[i * n_temp + j] = temp_list[j];

            status = flash_calculation_stability_analysis_QNSS(phase, K0, 1e-10);

            if (status == 1) {
                if (phase->phase_no == 0) {
                    sm->F[i * n_temp + j] = 0.0;
                }
                else {
                    sm->F[i * n_temp + j] = 1.0;
                }

                for (k = 0; k < ncomp; k++) {
                    sm->K[i * n_temp + j][k] = 0.0;
                }
            }
            else if (status == 0) {
                double *K00;

                K00 = malloc(ncomp * sizeof(*K00));

                if (fsa == NULL) {
                    Fv0 = Fv;

                    for(k = 0; k < ncomp; k++) {
                        if (K == NULL) {
                            K00[k] = K0[k];
                        }
                        else {
                            K00[k] = K[k];
                        }
                    }
                }
                else {
                    int flag = 0;
                    double input[ncomp + 2], K000[ncomp], prediction_time;

                    for (k = 0; k < ncomp; k++) {
                        input[k] = comp_X[k];
                    }
                    input[ncomp] = temp_list[j];
                    input[ncomp + 1] = pres_list[i];

#if 0
                    printf("input:");
                    for (k = 0; k < ncomp + 2; k++) {
                        printf("%e ", input[k]);
                    }
                    printf("\n");
#endif

                    prediction_time = flash_calculation_get_time(NULL);
                    flag = flash_calculation_split_ann_predict(fsa, input, ncomp + 2, 
                            &Fv0, K00);
                    prediction_time = flash_calculation_get_time(NULL) - prediction_time;
                    split_pred_time += prediction_time;

                    if (verb) {
                        printf("Prediction time: %e\n", prediction_time);
                    }

                    if (!flag) {
                        Fv0 = -1.0;
                    }
                    else {
                        if (Fv0 > 1.0) {
                            Fv0 = 0.9;
                        }

                        if (Fv0 < 0.0) {
                            Fv0 = 0.1;
                        }

#if 0
                        printf("Fv0: %e\n", Fv0);
                        printf("K: ");
#endif
                        flash_calculation_estimate_K(eos, K000);

                        for (k = 0; k < ncomp; k++) {
                            //printf("%e ", K00[k]);

                            if (K00[k] < 0.0 || fabs(log(K00[k])) < 1e-4) {
                                K00[k] = K000[k];
                            }
                        }
                        //printf("\n");
                    }
                }

                solve_time = flash_calculation_get_time(NULL);
                Fv = flash_calculation_two_phase_flash_Calculation_QNSS(eos, 
                        comp_X, K00, Fv0, 1e-10);
                solve_time = flash_calculation_get_time(NULL) - solve_time;
                split_solve_time += solve_time;

                if (verb && ((fabs(Fv) < 1e-5 || fabs(Fv - 1.0) < 1e-5))) {
                    printf("Split calculation time: %e\n", solve_time);
                }

#if 0
                printf("FFv: %e\n", Fv);
                printf("FK: \n");
#endif
                for(k = 0; k < ncomp; k++) {
                    K0[k] = K00[k];

#if 0
                    printf("%e ", K00[k]);
#endif
                }
                free(K00);
                //printf("\n");

                sm->F[i * n_temp + j] = Fv;
                for (k = 0; k < ncomp; k++) {
                    sm->K[i * n_temp + j][k] = K0[k];
                }

                K = sm->K[i * n_temp + j];
            }

            for (k = 0; k < ncomp; k++) {
                sm->x[i * n_temp + j][k] = comp_X[k];
            }
        }
    }

    if (output_name != NULL) {
        flash_calculation_output_split_calculation_map(sm, comp_X, 
                ncomp, 0, output_name);
    }

    free(pres_list);
    free(temp_list);
    flash_calculation_phase_free(&phase);
    free(eos);

    return sm;
}

SPLIT_PM_MAP * flash_calculation_draw_split_calculation_map_PM(COMP_LIST *comp_list, 
        double *comp_X, double T, double P_min, double P_max, double dP, 
        int selected_component, double dxx, FLASH_SPLIT_ANN *fsa, char *output_name)
{
    int i, j, k, ncomp = comp_list->ncomp;
    int n_pres, n_x;
    double *pres_list, *x_list;
    SPLIT_PM_MAP *sm;
    int status;
    double Fv, Fv0, *K0, *K, *x, sum_no_selected, 
           solve_time;
    EOS *eos;
    PHASE *phase;

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
    sm->n = n_pres * n_x;
    sm->pres = malloc(n_pres * n_x* sizeof(*(sm->pres)));
    sm->F = malloc(n_pres * n_x* sizeof(*(sm->F)));
    sm->K = malloc(n_pres * n_x* sizeof(*(sm->K)));
    for (i = 0; i < n_pres * n_x; i++) {
        sm->K[i] = malloc(ncomp * sizeof(*(sm->K[i])));
    }
    sm->x = malloc(n_pres * n_x* sizeof(*(sm->x)));
    for (i = 0; i < n_pres * n_x; i++) {
        sm->x[i] = malloc(ncomp * sizeof(*(sm->x[i])));
    }

    K0 = malloc(ncomp * sizeof(*K0));
    x = malloc(ncomp * sizeof(*x));
    eos = flash_calculation_EOS_new(comp_list, 0.0, T, 0);
    phase = flash_calculation_phase_new(eos, x);

    sum_no_selected = 0.0;
    for (i = 0; i < ncomp; i++) {
        if (i != selected_component) {
            sum_no_selected += comp_X[i];
        }
    }

    for (i = 0; i < n_pres; i++) {
        K = NULL;
        Fv = Fv0 = 0.0;

        for (j = 0; j < n_x; j++) {
            for (k = 0; k < ncomp; k++) {
                if (k == selected_component) {
                    x[k] = x_list[j];
                }
                else {
                    x[k] = (1.0 - x_list[j]) * comp_X[k] / sum_no_selected;
                }
            }

            for (k = 0; k < ncomp; k++) {
                sm->x[i * n_x + j][k] = x[k];
            }

            eos->pres = pres_list[i];
            eos->temp = T;

            sm->pres[i * n_x + j] = pres_list[i];

            status = flash_calculation_stability_analysis_QNSS(phase, K0, 1e-10);

            if (status == 1) {
                if (phase->phase_no == 0) {
                    sm->F[i * n_x + j] = 0.0;
                }
                else {
                    sm->F[i * n_x + j] = 1.0;
                }

                for (k = 0; k < ncomp; k++) {
                    sm->K[i * n_x + j][k] = 0.0;
                }
            }
            else if (status == 0) {
                double *K00;

                K00 = malloc(ncomp * sizeof(*K00));

                if (fsa == NULL) {
                    Fv0 = Fv;

                    for(k = 0; k < ncomp; k++) {
                        if (K == NULL) {
                            K00[k] = K0[k];
                        }
                        else {
                            K00[k] = K[k];
                        }
                    }
                }
                else {
                    int flag = 0;
                    double input[ncomp + 1], K000[ncomp], prediction_time;

                    for (k = 0; k < ncomp; k++) {
                        input[k] = x[k];
                    }
                    input[ncomp] = pres_list[i];

                    prediction_time = flash_calculation_get_time(NULL);
                    flag = flash_calculation_split_ann_predict(fsa, input, ncomp + 1, 
                            &Fv0, K00);
                    prediction_time = flash_calculation_get_time(NULL) - prediction_time;
                    split_pred_time += prediction_time;

                    if (verb) {
                        printf("Prediction time: %e\n", prediction_time);
                    }

                    if (!flag) {
                        Fv0 = -1.0;
                    }
                    else {
                        if (Fv0 > 1.0) {
                            Fv0 = 0.9;
                        }

                        if (Fv0 < 0.0) {
                            Fv0 = 0.1;
                        }

                        if (verb) {
                        printf("Fv0: %e\n", Fv0);
                        printf("K: ");
                        }

                        flash_calculation_estimate_K(eos, K000);

                        for (k = 0; k < ncomp; k++) {
                            if (verb) {
                            printf("%e ", K00[k]);
                            }

                            if (K00[k] < 0.0 || fabs(log(K00[k])) < 1e-4) {
                                K00[k] = K000[k];
                            }
                        }

                        if (verb) {
                        printf("\n");
                        }
                    }
                }

                solve_time = flash_calculation_get_time(NULL);
                Fv = flash_calculation_two_phase_flash_Calculation_QNSS(eos, 
                        x, K00, Fv0, 1e-10);
                solve_time = flash_calculation_get_time(NULL) - solve_time;
                split_solve_time += solve_time;

                if (verb && ((fabs(Fv) < 1e-5 || fabs(Fv - 1.0) < 1e-5))) {
                    printf("Split calculation time: %e\n", solve_time);
                }


                if (verb) {
                printf("Fv: %e\n", Fv);
                printf("KK: ");
                }
                for(k = 0; k < ncomp; k++) {
                    K0[k] = K00[k];

                    if (verb) {
                    printf("%e ", K00[k]);
                    }
                }
                free(K00);
                if (verb) {
                printf("\n");
                }

                sm->F[i * n_x + j] = Fv;
                for (k = 0; k < ncomp; k++) {
                    sm->K[i * n_x + j][k] = K0[k];
                }

                K = sm->K[i * n_x + j];
            }
        }
    }

    if (output_name != NULL) {
        flash_calculation_output_split_calculation_map_PM(sm, NULL, 
                ncomp, 0, output_name);
    }

    free(pres_list);
    free(x_list);
    free(K0);
    flash_calculation_phase_free(&phase);
    free(eos);

    return sm;
}

void flash_calculation_split_map_free(SPLIT_MAP **sm)
{
    int i;
    SPLIT_MAP *sm0 = *sm;

    free(sm0->temp);
    free(sm0->pres);
    free(sm0->F);

    for (i = 0; i < sm0->n; i++) {
        free(sm0->K[i]);
        free(sm0->x[i]);
    }

    free(sm0->K);

    free(*sm);
}

void flash_calculation_split_PM_map_free(SPLIT_PM_MAP **sm)
{
    int i;
    SPLIT_PM_MAP *sm0 = *sm;

    free(sm0->pres);
    free(sm0->F);

    for (i = 0; i < sm0->n; i++) {
        free(sm0->K[i]);
        free(sm0->x[i]);
    }

    free(sm0->K);
    free(sm0->x);

    free(*sm);
}

double flash_calculation_split_time_cost()
{
    return split_solve_time;
}

int flash_calculation_split_iteration_number()
{
    return split_itr;
}

int flash_calculation_split_failure_number()
{
    return split_failure;
}

double flash_calculation_split_pred_time_cost()
{
    return split_pred_time;
}



