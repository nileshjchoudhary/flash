#include "fc.h"
/* 
# ### Phase Diagram Construction
# The procedure of phase diagram construction is as follows:
# 1. staturation envelope construction
# 2. critical point calculation
# 3. calculate the temperature-pressure line at give $F_v$
*/

PHASE_DIAGRAM * flash_calculation_phase_diagram_construction(EOS *eos, double *z, 
        double P_start, double P_end, double T_start, double T_end, 
        double *Fv_list, int nF, PHASE_ENVELOPE *pe, CRITICAL_POINT *cp, 
        double dP, double dT)
{
    /* default value: dP 1.0, dT 1.0, *pe = NULL */
    /* The number of lines */
    int i, j, ith, ncomp = eos->ncomp, T_nstep, P_nstep, T_nstep2;
    PHASE_DIAGRAM *pd;
    double Tc, Pc;

    pd = malloc(sizeof(*pd));
    
    /* Calculate critical point */
    if (cp == NULL) {
        pd->cp = flash_calculation_critical_point_calculation(eos, z, 3.5, -1);
    }
    else {
        pd->cp = cp;
    }
    Tc = pd->cp->Tc;
    Pc = pd->cp->Pc;

    /* Construct phase envolope */
    if (pe == NULL) {
        pd->pe = flash_calculation_phase_saturation_envelope_construction(eos, 
                z, T_start, T_end, dT, 250.0, dP, P_end);
    }
    else {
        pd->pe = pe;
    }

    pd->pl = malloc(nF * sizeof(*(pd->pl)));
    pd->n_line = nF;
    
    T_nstep = pd->pe->n / 2;

    for (ith = 0; ith < nF; ith++) {
        double Fv, P_last, P, T;
        double *K, *K0;
        int first;
        double P_start0, T_start0;

        Fv = Fv_list[ith];

        pd->pl[ith] = malloc(sizeof(*(pd->pl[ith])));
        pd->pl[ith]->n = 0;
        pd->pl[ith]->F = Fv;
        pd->pl[ith]->P = malloc(T_nstep * sizeof(*(pd->pl[ith]->P)));
        pd->pl[ith]->T = malloc(T_nstep * sizeof(*(pd->pl[ith]->T)));

        P_last = 0.0;
        P = 0.0;
        K = NULL;
        K0 = malloc(ncomp * sizeof(*K0));
        
        /* Increase T and find the corresponding P */
        first = 1;
        for (i = 0; i < T_nstep; i++) { 
            T = pd->pe->Ts[i];

            if (P >= 1.0) {
                K = flash_calculation_calculate_status_at_fixed_temperature_and_Fv(eos, 
                        z, T, Fv, pd->pe->Ps[i], pd->pe->Ps[pd->pe->n - i - 1], P, K);
                P = eos->pres;
            }
            else {
                K = flash_calculation_calculate_status_at_fixed_temperature_and_Fv(eos, 
                        z, T, Fv, pd->pe->Ps[i], pd->pe->Ps[pd->pe->n - i - 1], -1, NULL);
                P = eos->pres;
            }
            
            if (P > 1.0) {
                int k;

                k = pd->pl[ith]->n;

                pd->pl[ith]->P[k] = P;
                pd->pl[ith]->T[k] = T;

                for (j = 0; j < ncomp; j++) {
                    K0[j] = K[j];
                }

                first = 0;

                pd->pl[ith]->n += 1;
            }
            
            if (fabs(P - Pc) < dP && fabs(T - Tc) < dT)
                break;
                
            if (T > Tc && !first)
                break;
        }

        if (K != NULL)  {
            free(K);
            K = NULL;
        }
                
        /* Check the distance between the last found point with critical
        # point. If the pressure distance or the temperature distance 
        # is larger than dP or dT, increase P and find the correspoinding 
        # T
        */
        P_nstep = 0;
        P_start0 = 0.;
        T_start0 = 0.;
        if (pd->pl[ith]->n == 0) {
            P_start0 = P_start;
            T_start0 = T_start;
        }
        else {
            P_start0 = pd->pl[ith]->P[pd->pl[ith]->n - 1];
            T_start0 = pd->pl[ith]->T[pd->pl[ith]->n - 1];
        }
            
        if (fabs(P_start0 - Pc) < dP)
            P_nstep = 0;
        else {
            P_nstep = (int)(fabs(P_start0 - Pc) / dP);
        }

        if (P_nstep > 0) {
            if (P_start0 > Pc)
                dP *= -1.0;
            
            T = 0.0;
            P = P_start0;
            K = NULL;
            
            for (i = 0; i < P_nstep; i++) {
                P += dP;

                if (T > 1.0) {
                    K = flash_calculation_calculate_status_at_fixed_pressure_and_Fv(eos, z, P, Fv,
                            1000.0, 1.0, T, K);
                    T = eos->temp;
                }
                else {
                    K = malloc(ncomp * sizeof(*K));
                    for (j = 0; j < ncomp; j++) {
                        K[j] = K0[j];
                    }

                    K = flash_calculation_calculate_status_at_fixed_pressure_and_Fv(eos, z, P, Fv,
                            1000.0, 1.0, T_start0, K);
                    T = eos->temp;
                }

                if (T > 1.0) {
                    int k;

                    k = pd->pl[ith]->n;

                    pd->pl[ith]->P[k] = P;
                    pd->pl[ith]->T[k] = T;

                    for (j = 0; j < ncomp; j++) {
                        K0[j] = K[j];
                    }

                    P_last = P;

                    pd->pl[ith]->n += 1;
                }
            }

            if (K != NULL) {
                free(K);
                K = NULL;
            }
        }


        T_nstep2 = 0;
        if (fabs(pd->pl[ith]->T[pd->pl[ith]->n - 1] - Tc) < dT) {
            T_nstep2 = 0;
        }
        else {
            T_nstep2 = (int)(fabs(pd->pl[ith]->T[pd->pl[ith]->n - 1] - Tc) / dT);
        }

        if (T_nstep2 > 0) {
            P = 0.0;
            T = pd->pl[ith]->T[pd->pl[ith]->n - 1];
            K = NULL;

            for (i = 0; i < T_nstep2; i++) {
                T -= dT;
                if (P > 1.0) {
                    K = flash_calculation_calculate_status_at_fixed_temperature_and_Fv(eos, z, T, Fv,
                            1000.0, 1.0, P, K);
                    P = eos->pres;
                }
                else {
                    K = malloc(ncomp * sizeof(*K));
                    for (j = 0; j < ncomp; j++) {
                        K[j] = K0[j];
                    }

                    K = flash_calculation_calculate_status_at_fixed_temperature_and_Fv(eos, z, T, Fv, 
                            1000.0, 1.0, (pd->pl[ith]->P[pd->pl[ith]->n - 1] + P_last) * 0.5, K);
                    P = eos->pres;
                }
                
                if (P > 1.0) {
                    int k;

                    k = pd->pl[ith]->n;

                    pd->pl[ith]->P[k] = P;
                    pd->pl[ith]->T[k] = T;

                    pd->pl[ith]->n += 1;
                }
            }

            if (K != NULL) {
                free(K);
                K = NULL;
            }
        }
        
        free(K0);
    }
    
    return pd;
}

/*
# ## Draw Phase Diagram
# For a feed composition, this function will draw a phase diagram at the given 
vapour mole fraction, temperature range and pressure range. For the example, 
the following figure shows the phase diagram at $F_v = [0.3, 0.5, 0.8]$ for 
the temperature in [200 K, 340 K] and the pressure in [20 atm, 80 atm].
*/

PHASE_DIAGRAM * flash_calculation_draw_phase_diagram(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double *Fv_list, int nF, double dT, double dP, PHASE_ENVELOPE *pe, 
        CRITICAL_POINT *cp, char *output_name)
{
    PHASE_DIAGRAM *pd;
    char file_name[100];
    FILE *fp;
    int i, j, k, ncomp = comp_list->ncomp;
    EOS *eos;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, eos_type);

    pd = flash_calculation_phase_diagram_construction(eos, comp_X, 
            P_min, P_max, T_min, T_max, Fv_list, nF, pe, cp, 
            dP, dT);

    flash_calculation_phase_envelope_output(pd->pe,
            comp_X, ncomp, output_name);

    if (output_name != NULL) {
        for (j = 0; j < pd->n_line; j++) {
            sprintf(file_name, "%s-phase-diagram-F%lf.csv", output_name, pd->pl[j]->F);
            fp = fopen(file_name, "a");

            for (i = 0; i < ncomp; i++) {
                fprintf(fp, "Component %d,", i + 1);
            }
            fprintf(fp, "Temperature,Pressure,Fv\n");

            for (i = 0; i < pd->pl[j]->n; i++) {
                for (k = 0; k < ncomp; k++) {
                    fprintf(fp, "%lf,", comp_X[k]);
                }

                fprintf(fp, "%lf,%lf,%lf\n", pd->pl[j]->T[i], pd->pl[j]->P[i], pd->pl[j]->F);
            }
            fclose(fp);
        }
    }

    printf("Critical Point: (%lf, %lf)\n", pd->cp->Tc, pd->cp->Pc);

    free(eos);

    return pd;
}

void flash_calculation_phase_diagram_free(PHASE_DIAGRAM **pd)
{
    int i;
    PHASE_DIAGRAM *pd0 = *pd;

    flash_calculation_phase_envelope_free(&(pd0->pe));
    flash_calculation_critical_point_free(&(pd0->cp));

    for (i = 0; i < pd0->n_line; i++) {
        free(pd0->pl[i]->P);
        free(pd0->pl[i]->T);
        free(pd0->pl[i]);
    }

    free(pd0->pl);

    free(*pd);
}
