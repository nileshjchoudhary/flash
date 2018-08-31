#include "fc.h"

static int flash_calculation_generate_simplex_number(double start, double end, 
            double dx, int ncomp)
{
    int i, sum, sub_sum;

    if (ncomp == 2) {
        sum = (int)((end - start) / dx + 0.5) + 1;

        return sum;
    }

    sub_sum = (int)((end - start) / dx + 0.5) + 1;
    //dx = (end - start) / sub_sum;

    sum = 0;
    for (i = 0; i < sub_sum; i++) {
        int sub_sum0;
        double x0;

        x0 = start + i * dx;
        sub_sum0 = flash_calculation_generate_simplex_number(start, 
                end - x0, dx, ncomp - 1);

        sum += sub_sum0;
    }

    return sum;
}


int flash_calculation_generate_simplex(double start, double end, 
        double dx, int ncomp, double ***x_list)
{
    int i, j, k, count, sum, sub_sum;
    double **x_list0;

    if (ncomp == 2) {
        sum = (int)((end - start) / dx + 0.5) + 1;
        //dx = (end - start) / sum;
        //printf("start: %e, end: %e, dx: %e\n", start, 
        //        end, dx);
        //printf("dx: %e, sum: %d\n", dx, sum);

        *x_list = malloc(sum * sizeof(**x_list));
        x_list0 = *x_list;

        for (i = 0; i < sum; i++) {
            x_list0[i] = malloc(2 * sizeof(double));

            x_list0[i][0] = start + i * dx;
            x_list0[i][1] = end - (start + i * dx);
        }

        return sum;
    }

    sum = flash_calculation_generate_simplex_number(start, end, 
            dx, ncomp);

    sub_sum = (int)((end - start) / dx + 0.5) + 1;
    //dx = (end - start) / sub_sum;
    //printf("DX: %e, sub_sum: %d\n", dx, sub_sum);

    *x_list = malloc(sum * sizeof(**x_list));
    x_list0 = *x_list;
    for (i = 0; i < sum; i++) {
        x_list0[i] = malloc(ncomp * sizeof(double));
    }

    count = 0;
    for (i = 0; i < sub_sum; i++) {
        int sub_sum0;
        double x0, **x_list_sub;

        x0 = start + i * dx;
        sub_sum0 = flash_calculation_generate_simplex(start, end - x0, 
                dx, ncomp - 1, &x_list_sub);

        for (j = 0; j < sub_sum0; j++) {
            x_list0[count][0] = x0;

            for (k = 0; k < ncomp - 1; k++) {
                x_list0[count][k + 1] = x_list_sub[j][k];
            }

            count++;
            free(x_list_sub[j]);
        }

        free(x_list_sub);
    }

    if (count != sum) {
        printf("Generate simplex data number is wrong!\n");
        exit(1);
    }

    return sum;
}

SET_NO_LIST * flash_calculation_generate_simplex_set_no(double start, 
        double end, double dx, int ncomp)
{
    int i, sum, sub_sum, count;
    SET_NO_LIST *set_no_list;

    if (ncomp == 2) {
        sum = (int)((end - start) / dx + 0.5) + 1;

        set_no_list = malloc(sizeof(*set_no_list));
        set_no_list->set_begin = 0;
        set_no_list->set_end = sum;
        set_no_list->set_init = 0;
        set_no_list->parent = NULL;
        set_no_list->next = NULL;
        set_no_list->previous = NULL;

        set_no_list->nchild = 0;
        set_no_list->child = NULL;

        return set_no_list;
    }

    sub_sum = (int)((end - start) / dx + 0.5) + 1;
    set_no_list = malloc(sizeof(*set_no_list));
    set_no_list->set_begin = 0;
    set_no_list->set_end = flash_calculation_generate_simplex_number(start, 
            end, dx, ncomp);
    set_no_list->nchild = sub_sum;
    set_no_list->child = malloc(sub_sum * sizeof(*(set_no_list->child)));

    for (i = 0; i < sub_sum; i++) {
        double x0;
        SET_NO_LIST *set_no_list0;

        x0 = start + i * dx;
        set_no_list0 = flash_calculation_generate_simplex_set_no(start, 
                end - x0, dx, ncomp - 1);

        set_no_list0->parent = set_no_list;
        set_no_list->child[i] = set_no_list0;
    }

    for (i = 0; i < sub_sum; i++) {
        if (i == 0) {
            set_no_list->child[i]->previous = NULL;
        }
        else {
            set_no_list->child[i]->previous = set_no_list->child[i - 1];
        }

        if (i != sub_sum - 1) {
            set_no_list->child[i]->next = set_no_list->child[i + 1];
        }
        else {
            set_no_list->child[i]->next = NULL;
        }
    }

    count = 0;
    for (i = 0; i < sub_sum; i++) {
        set_no_list->child[i]->set_begin += count;
        set_no_list->child[i]->set_end += count;

        count += set_no_list->child[i]->set_end
            - set_no_list->child[i]->set_begin;
    }

    if (set_no_list->parent == NULL) {
        set_no_list->set_init = 0;
    }

    for (i = 0; i < sub_sum; i++) {
        if (i == 0) {
            if (set_no_list->parent != NULL) {
                set_no_list->child[i]->set_init 
                    = set_no_list->parent->set_init;
            }
            else {
                set_no_list->child[i]->set_init 
                    = set_no_list->set_init;
            }
        }
        else {
            set_no_list->child[i]->set_init 
                = set_no_list->child[i-1]->set_begin;
        }
    }

    return set_no_list;
}

void flash_calculation_generate_simplex_set_no_free(SET_NO_LIST **set_no_list)
{
    int i;
    SET_NO_LIST *set_no_list0 = *set_no_list;

    if (set_no_list0->child != NULL) {
        for (i = 0; i != set_no_list0->nchild; i++) {
            flash_calculation_generate_simplex_set_no_free(&(set_no_list0->child[i])); 
        }
    }
    else {
        free(*set_no_list);
    }
}

void flash_calculation_generate_simplex_set_no_print(SET_NO_LIST *set_no_list)
{
    int i;

    printf("--------------------------\n");
    printf("set_begin: %d, set_end: %d, set_init: %d\n", set_no_list->set_begin,
            set_no_list->set_end, set_no_list->set_init);

    if (set_no_list->next != NULL) {
        printf("  next begin: %d, end: %d\n", set_no_list->next->set_begin,
                set_no_list->next->set_end);
    }

    if (set_no_list->previous != NULL) {
        printf("  previous begin: %d, end: %d\n", set_no_list->previous->set_begin,
                set_no_list->previous->set_end);
    }

    for (i = 0; i < set_no_list->nchild; i++) {
        printf("    children[%d]: \n", i);
        flash_calculation_generate_simplex_set_no_print(set_no_list->child[i]);
    }
}

static void flash_calculation_saturation_pressure_pre_order(SET_NO_LIST *set_no_list,
        EOS *eos, double **z, double T, double Ps_u_est, double Ps_l_est, double dP, 
        double P_max, double **xs, double *Ps_u, double *Ps_l)
{
    int i, k, ncomp = eos->ncomp;
    double P, P0, P1;

    if (set_no_list->child != NULL) {
        for (i = 0; i != set_no_list->nchild; i++) {
            Ps_u_est = Ps_u[set_no_list->child[i]->set_init];
            Ps_l_est = Ps_l[set_no_list->child[i]->set_init];
            flash_calculation_saturation_pressure_pre_order(set_no_list->child[i],
                    eos, z, T, Ps_u_est, Ps_l_est, dP, P_max, xs, Ps_u, Ps_l);
        }
    }
    else {
        /* upper */
        if (Ps_u_est > 1.0) {
            P = Ps_u_est;
        }
        else {
            P = 100.0;
        }

        P0 = P1 = P;
        for (i = set_no_list->set_begin; 
                i < set_no_list->set_end; i++) {
            xs[i] = malloc(ncomp * sizeof(*(xs[i])));
            //printf("-------------------upper\n");
            for (k = 0; k < ncomp; k++) {
                xs[i][k] = z[i][k];
            }

            if (P1 > P0) {
                P = flash_calculation_saturation_calculation(eos, xs[i],
                        T, P + 0.1, 0, dP, P_max);
            }
            else {
                P = flash_calculation_saturation_calculation(eos, xs[i],
                        T, P - 0.1, 0, dP, P_max);
            }

            P1 = P;
            P0 = P1;

            Ps_u[i] = P;
        }

        /* lower */
        P = Ps_l_est;
        for (i = set_no_list->set_begin; 
                i < set_no_list->set_end; i++) {
            //printf("-------------------lower\n");
            for (k = 0; k < ncomp; k++) {
                xs[i][k] = z[i][k];
                //printf("x[%d]: %e\n", k, xs[i][k]);
            }

            P = flash_calculation_saturation_calculation(eos, xs[i],
                    T, P - 0.1, 1, dP, P_max);

            Ps_l[i] = P;
        }

        for (i = set_no_list->set_begin;
                i < set_no_list->set_end; i++) {
            double temp;

            if (Ps_l[i] > Ps_u[i]) {
                temp = Ps_u[i];
                Ps_u[i] = Ps_l[i];
                Ps_l[i] = temp;
            }
        }
    }
}

PS_SIMPLEX_ISOTHERM * 
flash_calculation_saturation_pressure_simplex_isotherm(COMP_LIST *comp_list,
        double **z, int nz, SET_NO_LIST *set_no_list, double T, 
        double Ps_u_est, double Ps_l_est, double dP, double P_max)
{
    PS_SIMPLEX_ISOTHERM *ps;
    EOS *eos;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);

    ps = malloc(sizeof(*ps));

    ps->T = T;
    ps->n = nz;
    ps->xs = malloc(nz * sizeof(*(ps->xs)));
    ps->Ps_u = malloc(nz * sizeof(*(ps->Ps_u)));
    ps->Ps_l = malloc(nz * sizeof(*(ps->Ps_l)));

    ps->Ps_u[0] = Ps_u_est;
    ps->Ps_l[0] = Ps_l_est;

    flash_calculation_saturation_pressure_pre_order(set_no_list,
            eos, z, T, Ps_u_est, Ps_l_est, dP, P_max, ps->xs, 
            ps->Ps_u, ps->Ps_l);

    free(eos);

    return ps;
}

void
flash_calculation_saturation_pressure_simplex_isotherm_free(PS_SIMPLEX_ISOTHERM **ps)
{
    int i;
    PS_SIMPLEX_ISOTHERM *ps0 = *ps;

    free(ps0->Ps_u);
    free(ps0->Ps_l);
    for (i = 0; i < ps0->n; i++) {
        free(ps0->xs[i]);
    }
    free(ps0->xs);

    free(*ps);
}


void flash_calculation_saturation_pressure_simplex_isotherm_output(PS_SIMPLEX_ISOTHERM *ps, 
        int ncomp, char *output)
{
    int i, k;
    char file_name[100];
    FILE *fp;

    sprintf(file_name, "%s-simplex-PS-PM.csv", output);

    fp = fopen(file_name, "a");

    for (i = 0; i < ps->n; i++) {
        for (k = 0; k < ncomp; k++) {
            fprintf(fp, "%e,", ps->xs[i][k]);
        }

        fprintf(fp, "%e,", ps->Ps_u[i]);
        fprintf(fp, "%e\n", ps->Ps_l[i]);
    }
    fclose(fp);
}

PS_SIMPLEX_ISOTHERM *
flash_calculation_saturation_pressure_simplex_isotherm_data(COMP_LIST *comp_list,
    double T, double dx, double dP, double P_max, char *output)
{
    int nx;
    double **x_list;
    SET_NO_LIST *set_no_list;
    PS_SIMPLEX_ISOTHERM *ps_simplex;

    set_no_list = flash_calculation_generate_simplex_set_no(0.0, 
            1.0, dx, comp_list->ncomp);
    nx = flash_calculation_generate_simplex(0.0, 1.0, dx,
            comp_list->ncomp, &x_list);

    ps_simplex = flash_calculation_saturation_pressure_simplex_isotherm(comp_list,
            x_list, nx, set_no_list, T, 100.0, 1.0, dP, P_max);

    flash_calculation_saturation_pressure_simplex_isotherm_output(ps_simplex, 
            comp_list->ncomp, output);

    flash_calculation_generate_simplex_set_no_free(&set_no_list);

    return ps_simplex;
}


SPLIT_SIMPLEX_ISOTHERM *
flash_calculation_split_simplex_isotherm_data(COMP_LIST *comp_list,
    PS_SIMPLEX_ISOTHERM *ps, double dP, char *output)
{
    SPLIT_SIMPLEX_ISOTHERM *sp;

    sp = flash_calculation_split_simplex_isotherm(comp_list,
            ps, dP);
    flash_calculation_split_simplex_isotherm_output(sp, 
            comp_list->ncomp, output);

    return sp;
}

SPLIT_SIMPLEX_ISOTHERM *
flash_calculation_split_simplex_isotherm(COMP_LIST *comp_list,
        PS_SIMPLEX_ISOTHERM *ps, double dP)
{
    int i, j, k, final_n, ncomp = comp_list->ncomp;
    SPLIT_SIMPLEX_ISOTHERM *sp;
    EOS *eos;
    double *K0;
    PHASE *phase;

    sp = malloc(sizeof(*sp));
    sp->T = ps->T;

    final_n = 0;
    for (i = 0; i < ps->n; i++) {
        double Psu, Psl;

        Psu = ps->Ps_u[i];
        Psl = ps->Ps_l[i];

        if (Psu - Psl > 1e-3) {
            final_n++;
        }
    }

    sp->n = final_n;
    sp->nP = malloc(sp->n * sizeof(*(sp->nP)));
    sp->P = malloc(sp->n * sizeof(*(sp->P)));
    sp->Fv = malloc(sp->n * sizeof(*(sp->Fv)));
    sp->K = malloc(sp->n * sizeof(*(sp->K)));
    sp->xs = malloc(sp->n * sizeof(*(sp->xs)));

    final_n = 0;
    for (i = 0; i < ps->n; i++) {
        double Psu, Psl;
        int nP0;

        Psu = ps->Ps_u[i];
        Psl = ps->Ps_l[i];

        if (Psu - Psl < 1e-3)
            continue;

        nP0 = (int)((Psu - Psl) / dP);

        if (fabs(Psl + nP0 * dP - Psu) > 1e-3) {
            nP0++;
        }
        sp->P[final_n] = malloc(nP0 * sizeof(*(sp->P[final_n])));

        sp->nP[final_n] = nP0;
        for (j = 0; j < nP0; j++) {
            if (Psu - (Psl + j * dP) > 1e-3) {
                sp->P[final_n][j] = Psl + j * dP;
            }
            else {
                sp->P[final_n][j] = Psu;
            }
        }

        sp->Fv[final_n] = malloc(nP0 * sizeof(*(sp->Fv[final_n])));
        sp->K[final_n] = malloc(nP0 * sizeof(*(sp->K[final_n])));
        for (j = 0; j < nP0; j++) {
            sp->K[final_n][j] = malloc(ncomp 
                    * sizeof(*(sp->K[final_n][j])));
        }
        sp->xs[final_n] = malloc(ncomp 
                * sizeof(*(sp->xs[final_n])));
        for (j = 0; j < ncomp; j++) {
            sp->xs[final_n][j] = ps->xs[i][j];
        }

        final_n++;
    }

    K0 = malloc(ncomp * sizeof(*K0));
    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);

    for (i = 0; i < sp->n; i++) {
        double Fv = 0.5;

        for (j = 0; j < sp->nP[i]; j++) {
            eos->pres = sp->P[i][j];
            eos->temp = sp->T;

            if (j == 0) {
                phase = flash_calculation_phase_new(eos, sp->xs[i]);
                flash_calculation_stability_analysis_QNSS(phase, K0, 1e-10);

                Fv = flash_calculation_two_phase_flash_Calculation_QNSS(eos, 
                        sp->xs[i], K0, Fv, 1e-10);

                flash_calculation_phase_free(&phase);
            }
            else {
                Fv = flash_calculation_two_phase_flash_Calculation_QNSS(eos, 
                        sp->xs[i], K0, Fv, 1e-10);
            }

            sp->Fv[i][j] = Fv;
            for (k = 0; k < ncomp; k++) {
                sp->K[i][j][k] = K0[k];
            }
        }
    }

    free(eos);

    return sp;
}


void flash_calculation_split_simplex_isotherm_free(SPLIT_SIMPLEX_ISOTHERM **sp)
{
    int i, j;
    SPLIT_SIMPLEX_ISOTHERM *sp0 = *sp;

    for (i = 0; i < sp0->n; i++) {
        free(sp0->P[i]);
        free(sp0->Fv[i]);

        for (j = 0; j < sp0->nP[i]; j++) {
            free(sp0->K[i][j]);
        }
        free(sp0->xs[i]);
    }

    free(sp0->nP);

    free(*sp);
}

void flash_calculation_split_simplex_isotherm_output(SPLIT_SIMPLEX_ISOTHERM *sp, 
        int ncomp, char *output)
{
    int i, j, k;
    char file_name[100];
    FILE *fp;

    sprintf(file_name, "%s-simplex-SPLIT-PM.csv", output);

    fp = fopen(file_name, "a");

    for (i = 0; i < sp->n; i++) {
        for (j = 0; j < sp->nP[i]; j++) {
            for (k = 0; k < ncomp; k++) {
                fprintf(fp, "%e,", sp->xs[i][k]);
            }

            fprintf(fp, "%e,", sp->P[i][j]);
            fprintf(fp, "%e", sp->Fv[i][j]);

            for (k = 0; k < ncomp; k++) {
                fprintf(fp, ",%e", sp->K[i][j][k]);
            }

            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}

void flash_calculation_simplex_isotherm_data(COMP_LIST *comp_list, 
        double T, double dx, double dP, double P_max, char *output)
{
    PS_SIMPLEX_ISOTHERM *ps;
    SPLIT_SIMPLEX_ISOTHERM *sp;

    printf("==========================================\n");
    ps = flash_calculation_saturation_pressure_simplex_isotherm_data(comp_list,
            T, dx, dP, P_max, output);
    printf("Saturation pressure calculation is done!\n");
    printf("=========================================\n");
    printf("\n");

    printf("=========================================\n");
    sp = flash_calculation_split_simplex_isotherm_data(comp_list,
            ps, dP, output);
    printf("phase split calculation is done!\n");
    printf("=========================================\n");
    printf("\n");

    flash_calculation_saturation_pressure_simplex_isotherm_free(&ps);
    flash_calculation_split_simplex_isotherm_free(&sp);
}







