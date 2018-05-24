#include "fc.h"
/* ## 6. Phase Diagram Construction
# The following code is designed for phase diagram construntion.

# ### Calculate the values of equilibrium equations and Rachford-Rice equation
# The equilibrium equations and Rachford-Rice equation are the governing equations
*/
double flash_calculation_phase_diagram_construction_calculate_equations(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *G)
{
    int i, neqn = phase_L->ncomp + 1, ncomp = phase_L->ncomp;
	double error;
    
    /* Equilirium equations */
    flash_calculation_calculate_equilibrium_equation(phase_L, phase_V, G);
    
    /* Material equation */
    G[neqn - 1] = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_v, ncomp);
    
    error = 0.0;
    for (i = 0; i < neqn; i++) {
        error += G[i] * G[i];
	}
    
    return error;
}

/* ### Calculate the derivatives of equilibrium equations and Rachford-Rice equation
# The primary variables are K-values and pressure $P$.
*/
/*
        G[i] = ln(K[i]) + ln(phi_v[i](x_v)) - ln(phi_l[i](x_l))
        dG[i]/dK[j] = 1.0 / K[i] * sigma_ij
                    + 1.0 / phi_v[i](x_v) * dphi_v[i]/dx_v[j] * dx_v[j]/dK[j]
                    - 1.0 / phi_l[i](x_l) * dphi_l[i]/dx_l[i] * dx_l[j]/dK[j]
        dG[i]/dp = 1.0 / phi_v[i](x_v) * dphi_v[i]/dp
                 - 1.0 / phi_l[i](x_l) * dphi_l[i]/dp
        
        G[Nc+1] = Sum_i((K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0)))
        dG[Nc+1]/dK[i] = z[i] / (1.0 + F_v * (K[i] - 1.0))
                    - (K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0))^2 * F_v
        dG[Nc+1]/dp = 0.0
*/
void flash_calculation_phase_diagram_construction_calculate_equations_derivative(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *dG)
{
    int i, j, ncomp = phase_L->ncomp;
    double *dx_l, *dx_v;

    dx_l = malloc(ncomp * sizeof(*dx_l));
    dx_v = malloc(ncomp * sizeof(*dx_v));

    flash_calculation_calculate_composition_derivative(K, z, F_v, dx_l, dx_v, ncomp);

    /* dGi/dKj */
    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                dG[i * (ncomp + 1) + j] = 1.0 / K[i];
			}
            else {
                dG[i * (ncomp + 1) + j] = 0.0;
			}
            
            dG[i * (ncomp + 1) + j] += 1.0 / phase_V->phi[i] * phase_V->dphi_dx[i * ncomp + j] * dx_v[j]
				- 1.0 / phase_L->phi[i] * phase_L->dphi_dx[i * ncomp + j] * dx_l[j];
		}
	}

    /* dGi/dp */
    for (i = 0; i < ncomp; i++) {
		dG[i * (ncomp + 1) + ncomp] = 1.0 / phase_V->phi[i] * phase_V->dphi[i]
			- 1.0 / phase_L->phi[i] * phase_L->dphi[i];
	}
    
    /* dGnc+1/dKi */
    for (i = 0; i < ncomp; i++) {
		double temp;

        temp = 1.0 + F_v * (K[i] - 1.0);
        dG[ncomp * (ncomp + 1) + i] = z[i] / temp 
			- (K[i] - 1.0) * z[i] / (temp * temp) * F_v;
	}
    
    /* dGnc+1/dp */
    dG[ncomp * (ncomp + 1) + ncomp] = 0.0;

	free(dx_l);
	free(dx_v);
}

/*
        G[i] = ln(K[i]) + ln(phi_v[i](x_v)) - ln(phi_l[i](x_l))
        dG[i]/dK[j] = 1.0 / K[i] * sigma_ij
                    + 1.0 / phi_v[i](x_v) * dphi_v[i]/dx_v[j] * dx_v[j]/dK[j]
                    - 1.0 / phi_l[i](x_l) * dphi_l[i]/dx_l[i] * dx_l[j]/dK[j]
        dG[i]/dT = 1.0 / phi_v[i](x_v) * dphi_v[i]/dT
                 - 1.0 / phi_l[i](x_l) * dphi_l[i]/dT
        
        G[Nc+1] = Sum_i((K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0)))
        dG[Nc+1]/dK[i] = z[i] / (1.0 + F_v * (K[i] - 1.0))
                    - (K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0))^2 * F_v
        dG[Nc+1]/dp = 0.0
*/

void flash_calculation_phase_diagram_construction_calculate_equations_derivative_with_K_T(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *dG)
{
    int i, j, ncomp = phase_L->ncomp;
	double *dx_l, *dx_v;
    
    dx_l = malloc(ncomp * sizeof(*dx_l));
    dx_v = malloc(ncomp * sizeof(*dx_v));

    flash_calculation_calculate_composition_derivative(K, z, F_v, dx_l, dx_v, ncomp);

    /* dGi/dKj */
    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                dG[i * (ncomp + 1) + j] = 1.0 / K[i];
			}
            else {
                dG[i * (ncomp + 1) + j] = 0.0;
			}
            
            dG[i * (ncomp + 1) + j] += 1.0 / phase_V->phi[i] * phase_V->dphi_dx[i * ncomp + j] * dx_v[j]
				- 1.0 / phase_L->phi[i] * phase_L->dphi_dx[i * ncomp + j] * dx_l[j];
		}
	}
    
    /* dGi/dp */
    for (i = 0; i < ncomp; i++) {
        dG[i * (ncomp + 1) + ncomp] = 1.0 / phase_V->phi[i] * phase_V->dphi_dT[i]
			- 1.0 / phase_L->phi[i] * phase_L->dphi_dT[i];
	}
    
    /* dGnc+1/dKi */
    for (i = 0; i < ncomp; i++) {
		double temp;

        temp = 1.0 + F_v * (K[i] - 1.0);
        dG[ncomp * (ncomp + 1) + i] = z[i] / temp 
			- (K[i] - 1.0) * z[i] / (temp * temp) * F_v;
	}
    
    /* dGnc+1/dp */
    dG[ncomp * (ncomp + 1) + ncomp] = 0.0;

	free(dx_l);
	free(dx_v);
}

double * flash_calculation_phase_diagram_construction_update_variables(double *dG, double *G, int dim)
{
	int i;
	double *x;

	for (i = 0; i < dim; i++) {
		G[i] = -G[i];
	}

	x = malloc(dim * sizeof(*x));

	flash_calculation_solve_dense_linear_system(dG, G, x, dim);
    
    return x;
}

int flash_calculation_search_unstable_temperature(COMP_LIST *comp_list, 
        double *comp_X, double P, double *T_list, int nT)
{
    int i, status;
    EOS *eos;
    PHASE *phase;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);
    phase = flash_calculation_phase_new(eos, comp_X);

    for (i = 0; i < nT; i++) {
        eos->pres = P;
        eos->temp = T_list[i];

        status = flash_calculation_stability_analysis_QNSS(phase, NULL, 1e-10);
        
        if (status == 0) {
            break;
        }
    }

    flash_calculation_phase_free(&phase);
    free(eos);
            
    return i;
}

int flash_calculation_search_stable_temperature(COMP_LIST *comp_list, 
        double *comp_X, double P, double *T_list, int nT)
{
    int i, status;
    EOS *eos;
    PHASE *phase;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);
    phase = flash_calculation_phase_new(eos, comp_X);

    for (i = 0; i < nT; i++) {
        eos->pres = P;
        eos->temp = T_list[i];

        status = flash_calculation_stability_analysis_QNSS(phase, NULL, 1e-10);
        
        if (status == 1) {
            break;
        }
    }
            
    return i;
}


/* ### Saturation Envelope Construction
# We devide the temperature into some intervals wit $\delta T$. 
    At each $T$, we calculate the saturation pressure. 
    We first calculate the upper saturation pressure (bubble point pressure), 
    then the lower saturation pressure (dew point pressure).
*/
PHASE_ENVELOPE * flash_calculation_phase_saturation_envelope_construction(EOS *eos, 
        double *z, double T_start, double T_end, double dT, double P_est, double dP, 
        double P_max)
{
	double dT_min = 0.1, *T_list, P, last;
	int count, n, i, begin_index;
	int flag, flag2, flag3;
    PHASE_ENVELOPE *pe;

	n = (int)((T_end - T_start) / dT) + 1;
	T_list = malloc(n * sizeof(*T_list));

	for (i = 0; i < n - 1; i++) {
		T_list[i] = T_start + dT * i;
	}
	T_list[n - 1] = T_end;

    begin_index = flash_calculation_search_unstable_temperature(eos->comp_list, 
            z, 1.1, T_list, n);

    pe = malloc(sizeof(*pe));
    
    pe->Ps = malloc((2 * n + 20) * sizeof(*(pe->Ps)));
    pe->Ts = malloc((2 * n + 20) * sizeof(*(pe->Ts)));

	count = 0;
    flag = 0;
    P = P_est;

    for (i = 0; i < begin_index; i++) {
        pe->Ps[count] = 1.0;
        pe->Ts[count] = T_list[i];
        count += 1;
    }

    for (i = begin_index; i < n; i++) {
        //printf("----- Upper: temperature %lf\n", T_list[i]);
        P = flash_calculation_saturation_calculation(eos, z, T_list[i], 
                P, 0, dP, P_max);
        //printf("             pressure    %lf\n", P); 
        
        if (P > 1.0) {
            pe->Ps[count] = P;
            pe->Ts[count] = T_list[i];
            count += 1;
            flag = 1;
		}
        else {
            if (flag) {
				double T_l, T_r, T_c;

                T_l = T_list[i - 1];
                T_r = T_list[i];
                
                P = pe->Ps[count - 1];
                T_c = (T_l + T_r) * 0.5;

                while(1) {
                    P = flash_calculation_saturation_calculation(eos, z, T_c, P, 0, dP * 0.1,
                            P_max);
                    
                    if (P > 1.0) {
                        pe->Ps[count] = P;
                        pe->Ts[count] = T_c;
                        T_l = T_c;
                        count += 1;
					}
                    else {
                        T_r = T_c;
                        P = pe->Ps[count - 1];
					}
                        
                    if (fabs(T_r - T_l) < dT_min) {
                        break;
					}
                    else {
                        T_c = (T_l + T_r) * 0.5;
					}
				}

                break;
			}
            else {
                pe->Ps[count] = P;
                pe->Ts[count] = T_list[i];
                P = P_est;
                count += 1;
			}
		}
	}

    while (1) {
        int count0, *insert_index, insert_count, insert_count0;
        int flag_dP = 1;
        double *insert_P, *insert_T, *Ps_tmp, *Ts_tmp;

        count0 = count;

        insert_count = 0;
        insert_index = malloc(count0 * sizeof(*insert_index));
        insert_P = malloc(count0 * sizeof(*insert_P));
        insert_T = malloc(count0 * sizeof(*insert_T));

        for (i = 0; i < count0 - 1; i++) {
            if (pe->Ps[i] > 1.0 && pe->Ps[i + 1] > 1.0
                    && pe->Ps[i] < P_max) {
                if (fabs(pe->Ps[i] - pe->Ps[i + 1]) > dP
                        && fabs(pe->Ts[i] - pe->Ts[i + 1]) > dT * 0.01) {
                    flag_dP = 0;
                    
                    insert_index[insert_count] = i;

                    insert_T[insert_count] = (pe->Ts[i] + pe->Ts[i + 1]) * 0.5;
                    insert_P[insert_count] 
                        = flash_calculation_saturation_calculation(eos, z, 
                                insert_T[insert_count], pe->Ps[i], 0, dP * 0.1, 
                                P_max);

#if 0
                    printf("insert: Ps: %e %e, Ts: %e %e, insert_P: %e\n", pe->Ps[i], pe->Ps[i + 1],
                            pe->Ts[i], pe->Ts[i + 1], insert_P[insert_count]);
#endif

                    insert_count++;
                }
            }
        }

        //printf("insert_count: %d\n", insert_count);

        Ps_tmp = malloc(2 * (insert_count + count) * sizeof(double));
        Ts_tmp = malloc(2 * (insert_count + count) * sizeof(double));

        insert_count0 = 0;
        count0 = 0;
        for (i = 0; i < count; i++) {
            Ps_tmp[count0] = pe->Ps[i];
            Ts_tmp[count0] = pe->Ts[i];
            count0++;

            if (insert_count0 < insert_count 
                    && i == insert_index[insert_count0]) {
                Ps_tmp[count0] = insert_P[insert_count0];
                Ts_tmp[count0] = insert_T[insert_count0];
                count0++;
                insert_count0++;
            }
        }

        count += insert_count;

        free(pe->Ps);
        free(pe->Ts);

        pe->Ps = Ps_tmp;
        pe->Ts = Ts_tmp;

        free(insert_index);
        free(insert_P);
        free(insert_T);

        if (flag_dP) {
            break;
        }
    }


    flag = 0;
    flag2 = 0;
    flag3 = 0;
    P = pe->Ps[count - 1] - 1.0;
    if (P < 1.0)
        P = 1.0;
    last = count;

    for (i = last - 1; i > 0; i--) {
        //printf("----- Down : temperature %lf\n", pe->Ts[i]);
        if (flag2) {
            pe->Ps[count] = 1.0;
            pe->Ts[count] = pe->Ts[i];
			count += 1;
            continue;
		}

        if (flag3) {
            pe->Ps[count] = pe->Ps[i];
            pe->Ts[count] = pe->Ts[i];
			count += 1;
            continue;
		}
            
        P = flash_calculation_saturation_calculation(eos, z, pe->Ts[i], P, 1, dP * 0.1, P_max);
        //printf("             pressure    %lf\n", P);
        
        if (P > 1.0) {
            pe->Ps[count] = P;
            pe->Ts[count] = pe->Ts[i];
            flag = 1;
			count += 1;
		}
        else {
            if (flag) {
                pe->Ps[count] = 1.0;
                pe->Ts[count] = pe->Ts[i];
                flag2 = 1;
				count += 1;
			}
            else {
                P = pe->Ps[i - 1] - 5.0;
                if (P < 1.0)
                    P = 1.0;
            }
		}
	}

    pe->n = count;
    free(T_list);

    return pe;
}

void flash_calculation_phase_envelope_PM_output(PHASE_ENVELOPE_PM *pe_pm,
        int ncomp, int selected_component, char *output)
{
    char file_name[100], file_name_upper[100], file_name_down[100];
    int i, j, k;
    FILE *fp, *fp_upper, *fp_down;
    double xv, P;
    int flag;

    if (output == NULL)
        return;

    sprintf(file_name, "%s-envelope-PM.csv", output);
    fp = fopen(file_name, "a");

    for (k = 0; k < ncomp; k++) {
        fprintf(fp, "Component %d,", k);
    }
    fprintf(fp, "Pressure\n");

    for (i = 0; i < pe_pm->n; i++) {
        for (k = 0; k < ncomp; k++) {
            fprintf(fp, "%lf,", pe_pm->xs[i][k]);
        }

        fprintf(fp, "%lf\n", pe_pm->Ps[i]);
    }
    fclose(fp);

    sprintf(file_name_upper, "%s-envelope-PM-upper.csv", output);
    fp_upper = fopen(file_name_upper, "a");

    sprintf(file_name_down, "%s-envelope-PM-down.csv", output);
    fp_down = fopen(file_name_down, "a");

    flag = 0;
    for (i = 0; i < pe_pm->n; i++) {
        xv = pe_pm->xs[i][selected_component];
        P = pe_pm->Ps[i];

        if (i == 0 || (i > 0 && xv > pe_pm->xs[i - 1][selected_component])
                || (xv == pe_pm->xs[i - 1][selected_component] && flag)) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp_upper, "%lf,", pe_pm->xs[i][j]);
            }

            fprintf(fp_upper, "%lf\n", pe_pm->Ps[i]);

            flag = 0;
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp_down, "%lf,", pe_pm->xs[i][j]);
            }

            fprintf(fp_down, "%lf\n", pe_pm->Ps[i]);

            flag = 1;
        }
    }

    fclose(fp_upper);
    fclose(fp_down);
}

PHASE_ENVELOPE_PM * flash_calculation_phase_saturation_envelope_construction_PM(COMP_LIST *comp_list, 
        double *z, double T, double P_est, double dP, int selected_component, double dx, double P_max,
        char *output)
{
    int ncomp = comp_list->ncomp, i, k, n_x_list, count, last;
    double *x_list, *X, sum_no_selected, P;
    PHASE_ENVELOPE_PM *pe_pm;
    int flag, flag2, flag3 = 0;
    EOS *eos;
    FILE *fp;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);

    n_x_list = (int)((1.0 - 0.001) / dx);

    x_list = malloc(n_x_list * sizeof(*x_list));
    for (i = 0; i < n_x_list; i++) {
        x_list[i] = 0.001 + i * dx;
    }

    pe_pm = malloc(sizeof(*pe_pm));
    pe_pm->n = 0;
    pe_pm->Ps = malloc((2 * n_x_list + 20) * sizeof(*(pe_pm->Ps)));
    pe_pm->xs = malloc((2 * n_x_list + 20) * sizeof(*(pe_pm->xs)));

    sum_no_selected = 0.0;
    for (k = 0; k < ncomp; k++) {
        if (k != selected_component) {
            sum_no_selected += z[k];
        }
    }

    count = 0;
    flag = 0;
    P = P_est;

    X = malloc(ncomp * sizeof(*X));

    for (i = 0; i < n_x_list; i++) {
        //printf("X: ");
        for (k = 0; k < ncomp; k++) {
            if (k == selected_component) {
                X[k] = x_list[i];
            }
            else {
                X[k] = z[k] / sum_no_selected * (1.0 - x_list[i]);
            }

            //printf("%lf ", X[k]);
        }
#if 0
        printf("\n");
        printf("T: %lf\n", T);
        printf("Pest: %lf\n", P);
#endif

        P = flash_calculation_saturation_calculation(eos, X, T, P, 0, dP, P_max);
        //printf("Result: %lf\n", P);

        if (P > 1.0) {
            pe_pm->Ps[count] = P;
            pe_pm->xs[count] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                pe_pm->xs[count][k] = X[k];
            }
            count += 1;
            flag = 1;
        }
        else {
            if (flag) {
                double x_l, x_r, x_c;

                x_l = x_list[i - 1];
                x_r = x_list[i];

                P = pe_pm->Ps[count - 1];
                x_c = (x_l + x_r) * 0.5;

                for (k = 0; k < ncomp; k++) {
                    if (k == selected_component) {
                        X[k] = x_c;
                    }
                    else {
                        X[k] = z[k] / sum_no_selected * (1.0 - x_c);
                    }
                }

                while(1) {
                    P = flash_calculation_saturation_calculation(eos, X, T, P, 0, dP * 0.1, P_max);

                    if (P > 1.0) {
                        pe_pm->Ps[count] = P;
                        pe_pm->xs[count] = malloc(ncomp * sizeof(double));
                        for (k = 0; k < ncomp; k++) {
                            pe_pm->xs[count][k] = X[k];
                        }

                        x_l = x_c;
                        count += 1;
                    }
                    else {
                        x_r = x_c;
                        P = pe_pm->Ps[count - 1];
                    }

                    if (fabs(x_r - x_l) < 0.001) {
                        break;
                    }
                    else {
                        x_c = (x_l + x_r) * 0.5;

                        for (k = 0; k < ncomp; k++) {
                            if (k == selected_component) {
                                X[k] = x_c;
                            }
                            else {
                                X[k] = z[k] / sum_no_selected * (1.0 - x_c);
                            }
                        }
                    }
                }

                break;
            }
            else {
                pe_pm->Ps[count] = 1.0;
                pe_pm->xs[count] = malloc(ncomp * sizeof(double));
                for (k = 0; k < ncomp; k++) {
                    pe_pm->xs[count][k] = X[k];
                }
                count += 1;
            }
        }
    }

    free(X);

    flag = 0;
    flag2 = 0;
    if (pe_pm->Ps[count - 1] == 1.0) {
        flag2 = 1;
        flag3 = 1;
    }

    P = pe_pm->Ps[count - 1] - 1.0;
    if (P < 1.0)
        P = 1.0;
    last = count;

    for (i = last - 1; i >= 0; i--) {
        X = pe_pm->xs[i];

        //printf("----- Down : %d\n", i);
        if (flag2) {
            pe_pm->Ps[count] = 1.0;
            pe_pm->xs[count] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                pe_pm->xs[count][k] = X[k];
            }
            count += 1;
            continue;
        }

        P = flash_calculation_saturation_calculation(eos, X, T, P, 1, dP * 0.1, P_max);
        //printf("             pressure    %lf\n", P);

        if (P > 1.0) {
            pe_pm->Ps[count] = P;
            pe_pm->xs[count] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                pe_pm->xs[count][k] = X[k];
            }
            flag = 1;
            count += 1;
        }
        else {
            if (flag) {
                pe_pm->Ps[count] = 1.0;
                pe_pm->xs[count] = malloc(ncomp * sizeof(double));
                for (k = 0; k < ncomp; k++) {
                    pe_pm->xs[count][k] = X[k];
                }
                flag2 = 1;
                count += 1;
			}
            else {
                P = pe_pm->Ps[i - 1] - 5.0;
                if (P < 1.0) {
                    P = 1.0;
                }
            }
		}
    }

    pe_pm->n = count;

    if (flag3) {
        for (i = 0; i < count; i++) {
            free(pe_pm->xs[i]);
        }
        pe_pm->n = 0;
    }

    free(x_list);
    free(eos);

    flash_calculation_phase_envelope_PM_output(pe_pm,
            ncomp, selected_component, output);

    return pe_pm;
}

void flash_calculation_phase_envelope_output(PHASE_ENVELOPE *pe,
        double *comp_X, int ncomp, char *output_name)
{
    char file_name[100], file_name_upper[100], file_name_down[100];
    FILE *fp, *fp_upper, *fp_down;
    double T, P;
    int i, j;
    int flag, flag1;

    if (output_name == NULL) 
        return;

    sprintf(file_name, "%s-phase-envelope.csv", output_name);
    fp = fopen(file_name, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure\n");

    flag = 0; 
    flag1 = 0;
    for (i = 0; i < pe->n; i++) {
#if 0
        if (pe->Ps[i] > 1.0) {
            flag = 1;
        }

        if (!flag && !flag1 && i % 5 != 0) {
            flag1 = flag;

            continue;
        }
        flag1 = flag;
#endif

        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", pe->xs[i][j]);
            }
        }

        fprintf(fp, "%lf,%lf\n", pe->Ts[i], pe->Ps[i]);

    }
    fclose(fp);

    sprintf(file_name_upper, "%s-phase-envelope-upper.csv", output_name);
    fp_upper = fopen(file_name_upper, "a");

    sprintf(file_name_down, "%s-phase-envelope-down.csv", output_name);
    fp_down = fopen(file_name_down, "a");

    flag = 0; 
    flag1 = 0;
    for (i = 0; i < pe->n; i++) {
        T = pe->Ts[i];
        P = pe->Ps[i];

#if 0
        if (pe->Ps[i] > 1.0) {
            flag = 1;
        }
        else {
            flag = 0;
        }

        if (!flag && !flag1 && i % 5 != 0 
                && (i != pe->n - 1)) {
            flag1 = flag;

            continue;
        }
        flag1 = flag;

        if (i > 1) {
            if (fabs(T - pe->Ts[i - 1]) > 1e-5) {
                double dP, dT, rate, rate0;

                dP = P - pe->Ps[i - 1];
                dT = T - pe->Ts[i - 1];
                rate = fabs(dP / dT);

                dP = pe->Ps[i - 2] - pe->Ps[i - 1];
                dT = pe->Ts[i - 2] - pe->Ts[i - 1];

                if (dT == 0.0) {
                    rate0 = rate;
                }
                else {
                    rate0 = fabs(dP / dT);
                }

                if (rate > rate0 + 30.0 || rate + 30.0 < rate0) {
                    continue;
                }
            }
        }
#endif

        if (i == 0 || (i > 0 && T > pe->Ts[i - 1])) {
            if (comp_X != NULL) {
                for (j = 0; j < ncomp; j++) {
                    fprintf(fp_upper, "%lf,", comp_X[j]);
                }
            }
            else {
                for (j = 0; j < ncomp; j++) {
                    fprintf(fp_upper, "%lf,", pe->xs[i][j]);
                }
            }

            fprintf(fp_upper, "%lf,%lf\n", pe->Ts[i], pe->Ps[i]);
        }
        else {
            if (comp_X != NULL) {
                for (j = 0; j < ncomp; j++) {
                    fprintf(fp_down, "%lf,", comp_X[j]);
                }
            }
            else {
                for (j = 0; j < ncomp; j++) {
                    fprintf(fp_down, "%lf,", pe->xs[i][j]);
                }
            }

            fprintf(fp_down, "%lf,%lf\n", pe->Ts[i], pe->Ps[i]);
        }
    }

    fclose(fp_upper);
    fclose(fp_down);
}

void flash_calculation_phase_envelope_free(PHASE_ENVELOPE **pe)
{
    PHASE_ENVELOPE *pe0;

    pe0 = *pe;
    free(pe0->Ps);
    free(pe0->Ts);

    free(*pe);
}

void flash_calculation_phase_envelope_pm_free(PHASE_ENVELOPE_PM **pe_pm)
{
    int i;
    PHASE_ENVELOPE_PM *pe_pm0 = *pe_pm;

    free(pe_pm0->Ps);

    for (i = 0; i < pe_pm0->n; i++) {
        free(pe_pm0->xs[i]);
    }
    free(pe_pm0->xs);

    free(*pe_pm);
}
