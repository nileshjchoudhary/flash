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

static void flash_calculation_add_simplex_set_no(SET_NO_LIST *set_no_list,
        int count)
{
    int i;

    if (set_no_list->nchild == 0) {
        return;
    }

    for (i = 0; i < set_no_list->nchild; i++) {
        set_no_list->child[i]->set_begin += count;
        set_no_list->child[i]->set_end += count;

        flash_calculation_add_simplex_set_no(set_no_list->child[i], count);
    }

    return;
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
    set_no_list->parent = NULL;
    set_no_list->next = NULL;
    set_no_list->previous = NULL;
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

        if (i == sub_sum - 1) {
            set_no_list->child[i]->next = NULL;
        }
        else {
            set_no_list->child[i]->next = set_no_list->child[i + 1];
        }
    }

    count = 0;
    for (i = 0; i < sub_sum; i++) {
        set_no_list->child[i]->set_begin += count;
        set_no_list->child[i]->set_end += count;

        flash_calculation_add_simplex_set_no(set_no_list->child[i], count);

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
        EOS *eos, double **z, double *z_range, double T, double Ps_u_est, double Ps_l_est, 
        double dP, double P_max, double **xs, double *Ps_u, double *Ps_l)
{
    int i, k, ncomp = eos->ncomp;
    double P, P0, P1;

    if (set_no_list->child != NULL) {
        for (i = 0; i != set_no_list->nchild; i++) {
            Ps_u_est = Ps_u[set_no_list->child[i]->set_init];
            Ps_l_est = Ps_l[set_no_list->child[i]->set_init];
            flash_calculation_saturation_pressure_pre_order(set_no_list->child[i],
                    eos, z, z_range, T, Ps_u_est, Ps_l_est, dP, P_max, xs, Ps_u, Ps_l);
        }
    }
    else {
        int valid = 0;

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
            int flag = 1;

            xs[i] = malloc(ncomp * sizeof(*(xs[i])));

            for (k = 0; k < ncomp; k++) {
                if (z[i][k] < z_range[k*2] 
                        || z[i][k] > z_range[k*2+1]) {
                    flag = 0;
                    break;
                }
            }

            if (flag) {
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

                valid++;
            }
            else {
                for (k = 0; k < ncomp; k++) {
                    xs[i][k] = -1.0;
                }
                Ps_u[i] = -1.0;
            }
        }

        /* lower */
        P = Ps_l_est;
        for (i = set_no_list->set_begin; 
                i < set_no_list->set_end; i++) {
            int flag = 1;

            for (k = 0; k < ncomp; k++) {
                if (z[i][k] < z_range[k*2] 
                        || z[i][k] > z_range[k*2+1]) {
                    flag = 0;
                    break;
                }
            }

            if (flag) {
                //printf("-------------------lower\n");
                for (k = 0; k < ncomp; k++) {
                    xs[i][k] = z[i][k];
                    //printf("x[%d]: %e\n", k, xs[i][k]);
                }

                P = flash_calculation_saturation_calculation(eos, xs[i],
                        T, P - 0.1, 1, dP, P_max);

                Ps_l[i] = P;
            }
            else {
                for (k = 0; k < ncomp; k++) {
                    xs[i][k] = -1.0;
                }
                Ps_l[i] = -1.0;
            }
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

        if (valid > 0) {
            printf("    From %d to %d ... ", 
                    set_no_list->set_begin,
                    set_no_list->set_end);
            printf(" %d/%d are calculated.\n", valid, 
                    set_no_list->set_end - set_no_list->set_begin);
        }
    }
}

static void flash_calculation_saturation_pressure_bisection(double Pu_l, double Pu_r, 
        double Pl_l, double Pl_r, double *xs_l, double *xs_r,
        int set_no1, int set_no2,
        EOS *eos, double T, double dP, double P_max, PS_SIMPLEX_ISOTHERM *ps)
{
    int i, ncomp = eos->ncomp;

    if (Pu_l > 1.0 && Pu_r > 1.0 
            && (fabs(Pu_l - Pu_r) > dP 
            || fabs(Pl_l - Pl_r) > dP)) {
        double Pu_est, Pl_est, P;
        double Pu_m, Pl_m, *xs_m;

        Pu_est = (Pu_l + Pu_r) * 0.5;
        Pl_est = (Pl_l + Pl_r) * 0.5;

        ps->set_no[ps->n_adv] = malloc(2 * sizeof(int));
        ps->set_no[ps->n_adv][0] = set_no1;
        ps->set_no[ps->n_adv][1] = set_no2;

        ps->xs_adv[ps->n_adv] = malloc(ncomp * sizeof(double));
        for (i = 0; i < ncomp; i++) {
            ps->xs_adv[ps->n_adv][i] = (xs_l[i] + xs_r[i]) * 0.5;
        }
        xs_m = ps->xs_adv[ps->n_adv];


        P = flash_calculation_saturation_calculation(eos, 
                ps->xs_adv[ps->n_adv], T, Pu_est, 0,
                dP, P_max);
        ps->Ps_u_adv[ps->n_adv] = P;
        Pu_m = P;

        P = flash_calculation_saturation_calculation(eos, 
                ps->xs_adv[ps->n_adv], T, Pl_est, 1,
                dP, P_max);
        ps->Ps_l_adv[ps->n_adv] = P;
        Pl_m = P;

        ps->n_adv++;

        if (ps->n_adv > ps->alloc_adv - 1) {
            ps->alloc_adv += 20;
            ps->Ps_u_adv = realloc(ps->Ps_u_adv, 
                    ps->alloc_adv * sizeof(double));
            ps->Ps_l_adv = realloc(ps->Ps_l_adv, 
                    ps->alloc_adv * sizeof(double));
            ps->xs_adv = realloc(ps->xs_adv,
                    ps->alloc_adv * sizeof(*(ps->xs_adv)));
            ps->set_no = realloc(ps->set_no,
                    ps->alloc_adv * sizeof(*(ps->set_no)));
        }

        flash_calculation_saturation_pressure_bisection(Pu_l, Pu_m, 
                Pl_l, Pl_m, xs_l, xs_m, set_no1, set_no2, 
                eos, T, dP, P_max, ps);

        flash_calculation_saturation_pressure_bisection(Pu_m, Pu_r, 
                Pl_m, Pl_r, xs_m, xs_r, set_no1, set_no2,
                eos, T, dP, P_max, ps);
    }
}

static void flash_calculation_saturation_pressure_adaptive(SET_NO_LIST *set_no_list,
        EOS *eos, double T, double dP, double P_max, PS_SIMPLEX_ISOTHERM *ps)
{
    int i;
    SET_NO_LIST *next;

    if (set_no_list->parent == NULL) {
        ps->alloc_adv = 100;
        ps->Ps_u_adv = malloc(ps->alloc_adv * sizeof(double));
        ps->Ps_l_adv = malloc(ps->alloc_adv * sizeof(double));
        ps->xs_adv = malloc(ps->alloc_adv * sizeof(*(ps->xs_adv)));
        ps->set_no = malloc(ps->alloc_adv * sizeof(*(ps->set_no)));

        ps->n_adv = 0;
    }

    /* compare with its next neighbour */
    if (set_no_list->next != NULL) {
        double Pu, Pu_next, Pl, Pl_next;
        double *xs, *xs_next;

        next = set_no_list->next;

        Pu = ps->Ps_u[set_no_list->set_begin];
        Pu_next = ps->Ps_u[next->set_begin];

        Pl = ps->Ps_l[set_no_list->set_begin];
        Pl_next = ps->Ps_l[next->set_begin];

        xs = ps->xs[set_no_list->set_begin];
        xs_next = ps->xs[next->set_begin];

        flash_calculation_saturation_pressure_bisection(Pu, Pu_next, 
                Pl, Pl_next, xs, xs_next, set_no_list->set_begin,
                next->set_begin, eos, T, dP, P_max, ps);
    }

    for (i = set_no_list->set_begin; i < set_no_list->set_end - 1; i++) {
        double Pu, Pu_next, Pl, Pl_next;
        double *xs, *xs_next;

        Pu = ps->Ps_u[i];
        Pu_next = ps->Ps_u[i + 1];

        Pl = ps->Ps_l[i];
        Pl_next = ps->Ps_l[i + 1];

        xs = ps->xs[i];
        xs_next = ps->xs[i + 1];

        flash_calculation_saturation_pressure_bisection(Pu, Pu_next, 
                Pl, Pl_next, xs, xs_next, i, i + 1, eos, T, dP, P_max, ps);
    }



    /* go to its children */
    if (set_no_list->child != NULL) {
        for (i = 0; i < set_no_list->nchild; i++) {
            flash_calculation_saturation_pressure_adaptive(set_no_list->child[i],
                    eos, T, dP, P_max, ps);
        }
    }
}

    PS_SIMPLEX_ISOTHERM * 
flash_calculation_saturation_pressure_simplex_isotherm(COMP_LIST *comp_list,
        double **z, int nz, double *z_range, SET_NO_LIST *set_no_list, double T, 
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
            eos, z, z_range, T, Ps_u_est, Ps_l_est, dP, P_max, ps->xs, 
            ps->Ps_u, ps->Ps_l);

    printf("Adaptively adding points: \n");
    flash_calculation_saturation_pressure_adaptive(set_no_list,
            eos, T, dP, P_max, ps);

    printf("Adaptive: %d points are added.\n", ps->n_adv);

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

    free(ps0->Ps_u_adv);
    free(ps0->Ps_l_adv);
    for (i = 0; i < ps0->n_adv; i++) {
        free(ps0->xs_adv[i]);
        free(ps0->set_no[i]);
    }
    free(ps0->xs_adv);
    free(ps0->set_no);

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
        int flag = 1;

        for (k = 0; k < ncomp; k++) {
            if (ps->xs[i][k] < 0.0) {
                flag = 0;
                break;
            }
        }

        if (flag) {
            if (ps->Ps_u[i] > 1.0 || ps->Ps_l[i] > 1.0) {
                for (k = 0; k < ncomp; k++) {
                    fprintf(fp, "%e,", ps->xs[i][k]);
                }

                fprintf(fp, "%e,", ps->Ps_u[i]);
                fprintf(fp, "%e\n", ps->Ps_l[i]);
            }
        }
    }

    for (i = 0; i < ps->n_adv; i++) {
        if (ps->Ps_u_adv[i] > 1.0 || ps->Ps_l_adv[i] > 1.0) {
            for (k = 0; k < ncomp; k++) {
                fprintf(fp, "%e,", ps->xs_adv[i][k]);
            }

            fprintf(fp, "%e,", ps->Ps_u_adv[i]);
            fprintf(fp, "%e\n", ps->Ps_l_adv[i]);
        }
    }

    fclose(fp);
}

PS_SIMPLEX_ISOTHERM *
flash_calculation_saturation_pressure_simplex_isotherm_data(double **x_list, int nx,
        SET_NO_LIST *set_no_list, COMP_LIST *comp_list,
        double T, double *z_range, double dP, double P_max, char *output)
{
    PS_SIMPLEX_ISOTHERM *ps_simplex;

    ps_simplex = flash_calculation_saturation_pressure_simplex_isotherm(comp_list,
            x_list, nx, z_range, set_no_list, T, 100.0, 1.0, dP, P_max);

    flash_calculation_saturation_pressure_simplex_isotherm_output(ps_simplex, 
            comp_list->ncomp, output);


    return ps_simplex;
}


SPLIT_SIMPLEX_ISOTHERM *
flash_calculation_split_simplex_isotherm_data(SET_NO_LIST *set_no_list,
        COMP_LIST *comp_list,
        PS_SIMPLEX_ISOTHERM *ps, double dP, double P_min, double P_max,
        char *output)
{
    SPLIT_SIMPLEX_ISOTHERM *sp;

    sp = flash_calculation_split_simplex_isotherm(set_no_list, comp_list,
            ps, dP, P_min, P_max);
    flash_calculation_split_simplex_isotherm_output(sp, 
            comp_list->ncomp, output);

    return sp;
}

static void flash_calculation_split_simplex_pre_order(SET_NO_LIST *set_no_list,
        EOS *eos, int *nP, double **P, double **Fv, double ***K, double **xs, 
        double T, double Fv_est, double *K_est) 
{
    int i, j, k, ncomp = eos->ncomp;

    if (set_no_list->child != NULL) {
        for (i = 0; i < set_no_list->nchild; i++) {
            if (nP[set_no_list->child[i]->set_init] > 0) {
                Fv_est = Fv[set_no_list->child[i]->set_init][0];
                K_est = K[set_no_list->child[i]->set_init][0];
                flash_calculation_split_simplex_pre_order(set_no_list->child[i],
                        eos, nP, P, Fv, K, xs, T, Fv_est, K_est);
            }
            else {
                Fv_est = -1.0;
                K_est = NULL;
                flash_calculation_split_simplex_pre_order(set_no_list->child[i],
                        eos, nP, P, Fv, K, xs, T, Fv_est, K_est);
            }
        }
    }
    else {
        double *K0, Fv0 = 0.5;

        K0 = malloc(ncomp * sizeof(*K0));

        for (i = set_no_list->set_begin; 
                i < set_no_list->set_end; i++) {
            for (j = 0; j < nP[i]; j++) {
                eos->pres = P[i][j];
                eos->temp = T;

                if (j == 0 && i == set_no_list->set_begin) {
                    if (Fv_est <= 0.0 || Fv_est >= 1.0) {
                        PHASE *phase;

                        Fv0 = 0.5;

                        phase = flash_calculation_phase_new(eos, xs[i]);
                        flash_calculation_stability_analysis_QNSS(phase, K0, 1e-10);

                        flash_calculation_phase_free(&phase);
                    }
                    else {
                        for (k = 0; k < ncomp; k++) {
                            K0[k] = K_est[k];
                        }
                    }
                }

                if (Fv0 <= 0.0 || Fv0 >= 1.0) {
                    PHASE *phase;

                    if (Fv0 <= 0.0) 
                        Fv0 = 0.1;

                    if (Fv0 >= 1.0) 
                        Fv0 = 0.9;

                    phase = flash_calculation_phase_new(eos, xs[i]);
                    flash_calculation_stability_analysis_QNSS(phase, K0, 1e-10);

                    flash_calculation_phase_free(&phase);
                }

                Fv0 = flash_calculation_two_phase_flash_Calculation_QNSS(eos, 
                        xs[i], K0, Fv0, 1e-10);

                Fv[i][j] = Fv0;
                for (k = 0; k < ncomp; k++) {
                    K[i][j][k] = K0[k];
                }
            }
        }

        free(K0);
    }
}

SPLIT_SIMPLEX_ISOTHERM *
flash_calculation_split_simplex_isotherm(SET_NO_LIST *set_no_list, 
        COMP_LIST *comp_list,
        PS_SIMPLEX_ISOTHERM *ps, double dP, double P_min, double P_max)
{
    int i, j, k, ncomp = comp_list->ncomp;
    SPLIT_SIMPLEX_ISOTHERM *sp;
    EOS *eos;
    double *K0;
    PHASE *phase;

    sp = malloc(sizeof(*sp));
    sp->T = ps->T;

    sp->n = ps->n;
    sp->nP = malloc(sp->n * sizeof(*(sp->nP)));
    sp->P = malloc(sp->n * sizeof(*(sp->P)));
    sp->Fv = malloc(sp->n * sizeof(*(sp->Fv)));
    sp->K = malloc(sp->n * sizeof(*(sp->K)));
    sp->xs = malloc(sp->n * sizeof(*(sp->xs)));

    for (i = 0; i < ps->n; i++) {
        double Psu, Psl;
        int nP0;

        Psu = ps->Ps_u[i];
        Psl = ps->Ps_l[i];

        if (Psu - Psl < 1e-3 || Psu <= 0.0 || Psl <= 0.0) {
            sp->nP[i] = 0;

            continue;
        }

        if (Psu > P_max) {
            Psu = P_max;
        }
        if (Psl < P_min) {
            Psl = P_min;
        }

        nP0 = (int)((Psu - Psl) / dP);

        if (fabs(Psl + nP0 * dP - Psu) > 1e-3) {
            nP0++;
        }
        sp->P[i] = malloc(nP0 * sizeof(*(sp->P[i])));

        sp->nP[i] = nP0;
        for (j = 0; j < nP0; j++) {
            if (Psu - (Psl + j * dP) > 1e-3) {
                sp->P[i][j] = Psl + j * dP;
            }
            else {
                sp->P[i][j] = Psu;
            }
        }

        sp->Fv[i] = malloc(nP0 * sizeof(*(sp->Fv[i])));
        sp->K[i] = malloc(nP0 * sizeof(*(sp->K[i])));
        for (j = 0; j < nP0; j++) {
            sp->K[i][j] = malloc(ncomp 
                    * sizeof(*(sp->K[i][j])));
        }
        sp->xs[i] = malloc(ncomp 
                * sizeof(*(sp->xs[i])));
        for (j = 0; j < ncomp; j++) {
            sp->xs[i][j] = ps->xs[i][j];
        }
    }

    K0 = malloc(ncomp * sizeof(*K0));
    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);


    flash_calculation_split_simplex_pre_order(set_no_list,
            eos, sp->nP, sp->P, sp->Fv, sp->K, sp->xs, 
            sp->T, -1.0, NULL);

    sp->n_adv = ps->n_adv;
    sp->nP_adv = malloc(sp->n_adv * sizeof(*(sp->nP_adv)));
    sp->P_adv = malloc(sp->n_adv * sizeof(*(sp->P_adv)));
    sp->Fv_adv = malloc(sp->n_adv * sizeof(*(sp->Fv_adv)));
    sp->K_adv = malloc(sp->n_adv * sizeof(*(sp->K_adv)));
    sp->xs_adv = malloc(sp->n_adv * sizeof(*(sp->xs_adv)));

    for (i = 0; i < ps->n_adv; i++) {
        double Psu, Psl;
        int nP0;

        Psu = ps->Ps_u_adv[i];
        Psl = ps->Ps_l_adv[i];

        if (Psu - Psl < 1e-3 || Psu <= 0.0 || Psl <= 0.0)
            continue;

        if (Psu > P_max) {
            Psu = P_max;
        }
        if (Psl < P_min) {
            Psl = P_min;
        }

        nP0 = (int)((Psu - Psl) / dP);

        if (fabs(Psl + nP0 * dP - Psu) > 1e-3) {
            nP0++;
        }
        sp->P_adv[i] = malloc(nP0 * sizeof(*(sp->P_adv[i])));

        sp->nP_adv[i] = nP0;
        for (j = 0; j < nP0; j++) {
            if (Psu - (Psl + j * dP) > 1e-3) {
                sp->P_adv[i][j] = Psl + j * dP;
            }
            else {
                sp->P_adv[i][j] = Psu;
            }
        }

        sp->Fv_adv[i] = malloc(nP0 * sizeof(*(sp->Fv_adv[i])));
        sp->K_adv[i] = malloc(nP0 * sizeof(*(sp->K_adv[i])));
        for (j = 0; j < nP0; j++) {
            sp->K_adv[i][j] = malloc(ncomp 
                    * sizeof(*(sp->K_adv[i][j])));
        }
        sp->xs_adv[i] = malloc(ncomp * sizeof(*(sp->xs_adv[i])));
        for (j = 0; j < ncomp; j++) {
            sp->xs_adv[i][j] = ps->xs_adv[i][j];
        }
    }

    for (i = 0; i < sp->n_adv; i++) {
        double Fv = -1.0;

        for (j = 0; j < sp->nP_adv[i]; j++) {
            eos->pres = sp->P_adv[i][j];
            eos->temp = sp->T;

            if (j == 0) {
                int set_no1, set_no2;

                set_no1 = ps->set_no[i][0];
                set_no2 = ps->set_no[i][1];

                Fv = (sp->Fv[set_no1][0] + sp->Fv[set_no2][0]) * 0.5;

                for (k = 0; k < ncomp; k++) {
                    K0[k] = (sp->K[set_no1][0][k] 
                            + sp->K[set_no2][0][k]) * 0.5;
                }
            }

            if (Fv <= 0.0 || Fv >= 1.0) {
                if (Fv <= 0.0) 
                    Fv = 0.1;

                if (Fv >= 1.0) 
                    Fv = 0.9;

                phase = flash_calculation_phase_new(eos, sp->xs_adv[i]);
                flash_calculation_stability_analysis_QNSS(phase, K0, 1e-10);

                flash_calculation_phase_free(&phase);
            }

            Fv = flash_calculation_two_phase_flash_Calculation_QNSS(eos, 
                    sp->xs_adv[i], K0, Fv, 1e-10);

            sp->Fv_adv[i][j] = Fv;
            for (k = 0; k < ncomp; k++) {
                sp->K_adv[i][j][k] = K0[k];
            }
        }
    }

    free(K0);
    free(eos);

    return sp;
}


void flash_calculation_split_simplex_isotherm_free(SPLIT_SIMPLEX_ISOTHERM **sp)
{
    int i, j;
    SPLIT_SIMPLEX_ISOTHERM *sp0 = *sp;

    for (i = 0; i < sp0->n; i++) {
        if (sp0->nP[i] > 0) {
            free(sp0->P[i]);
            free(sp0->Fv[i]);

            for (j = 0; j < sp0->nP[i]; j++) {
                free(sp0->K[i][j]);
            }
            free(sp0->xs[i]);
        }
    }

    free(sp0->nP);

    for (i = 0; i < sp0->n_adv; i++) {
        free(sp0->P_adv[i]);
        free(sp0->Fv_adv[i]);

        for (j = 0; j < sp0->nP_adv[i]; j++) {
            free(sp0->K_adv[i][j]);
        }
        free(sp0->xs_adv[i]);
    }

    free(sp0->nP_adv);

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
            if (sp->Fv[i][j] > 1e-15 && sp->Fv[i][j] < 1.0 - 1e-15) {
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
    }

    for (i = 0; i < sp->n_adv; i++) {
        for (j = 0; j < sp->nP_adv[i]; j++) {
            if (sp->Fv_adv[i][j] > 1e-15 && sp->Fv_adv[i][j] < 1.0 - 1e-15) {
                for (k = 0; k < ncomp; k++) {
                    fprintf(fp, "%e,", sp->xs_adv[i][k]);
                }

                fprintf(fp, "%e,", sp->P_adv[i][j]);
                fprintf(fp, "%e", sp->Fv_adv[i][j]);

                for (k = 0; k < ncomp; k++) {
                    fprintf(fp, ",%e", sp->K_adv[i][j][k]);
                }

                fprintf(fp, "\n");
            }
        }
    }

    fclose(fp);
}

void flash_calculation_simplex_isotherm_data(COMP_LIST *comp_list, 
        double T, double dx, double *z_range, double dP, double P_min, 
        double P_max, char *output)
{
    int nx;
    double **x_list;
    SET_NO_LIST *set_no_list;
    PS_SIMPLEX_ISOTHERM *ps;
    SPLIT_SIMPLEX_ISOTHERM *sp;

    set_no_list = flash_calculation_generate_simplex_set_no(0.0, 
            1.0, dx, comp_list->ncomp);
    nx = flash_calculation_generate_simplex(0.0, 1.0, dx,
            comp_list->ncomp, &x_list);

    printf("Total points: %d\n", nx);

    printf("==========================================\n");
    ps = flash_calculation_saturation_pressure_simplex_isotherm_data(
            x_list, nx, set_no_list, comp_list,
            T, z_range, dP, P_max, output);
    printf("Saturation pressure calculation is done!\n");
    printf("=========================================\n");
    printf("\n");

    printf("=========================================\n");
    sp = flash_calculation_split_simplex_isotherm_data(
            set_no_list, comp_list,
            ps, dP, P_min, P_max, output);
    printf("phase split calculation is done!\n");
    printf("=========================================\n");
    printf("\n");

    flash_calculation_generate_simplex_set_no_free(&set_no_list);

    flash_calculation_saturation_pressure_simplex_isotherm_free(&ps);
    flash_calculation_split_simplex_isotherm_free(&sp);
}







