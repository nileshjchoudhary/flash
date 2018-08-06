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
        printf("start: %e, end: %e, dx: %e\n", start, 
                end, dx);
        printf("dx: %e, sum: %d\n", dx, sum);

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
    printf("DX: %e, sub_sum: %d\n", dx, sub_sum);

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
    double P;

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
        P = Ps_u_est;
        for (i = set_no_list->set_begin; 
                i < set_no_list->set_end; i++) {
            xs[i] = malloc(ncomp * sizeof(*(xs[i])));
            for (k = 0; k < ncomp; k++) {
                xs[i][k] = z[i][k];
                printf("xs[%d]: %e\n", k, xs[i][k]);
            }

            P = flash_calculation_saturation_calculation(eos, xs[i],
                    T, P, 0, dP, P_max);

            printf("upper saturation: %e\n", P);

            Ps_u[i] = P;
        }

        /* lower */
        P = Ps_l_est;
        for (i = set_no_list->set_begin; 
                i < set_no_list->set_end; i++) {
            for (k = 0; k < ncomp; k++) {
                xs[i][k] = z[i][k];
            }

            P = flash_calculation_saturation_calculation(eos, xs[i],
                    T, P, 1, dP, P_max);
            printf("lower saturation: %e\n", P);

            Ps_l[i] = P;
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

