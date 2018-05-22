#include "fc.h"

static int flash_calculation_generate_x_number(double start,
        double end, double dx, int ncomp)
{
    int i, sum, sum2;

    if (ncomp == 2) {
        sum = (int)((end - start) / dx);

        return sum + 1;
    }

    sum2 = (int)((end - start) / dx);
    sum = 0;
    for (i = 0; i < sum2; i++) {
        int sum3;
        double x0;

        x0 = start + i * dx;
        sum3 = flash_calculation_generate_x_number(start, 
                end - x0, dx, ncomp - 1);

        sum += sum3;
    }

    sum += ncomp - 1;

    return sum;
}

static int flash_calculation_generate_x(double start, double end, double dx, 
        int ncomp, double min_x, double ***x_list)
{
    int i, j, k, count, sum, sum2, ncomp2;
    double **x_list0;

    if (ncomp == 2) {
        sum = (int)((end - start) / dx);

        *x_list = malloc((sum + 1)* sizeof(**x_list));
        x_list0 = *x_list;

        //printf("N222222: %lf, %lf\n", start, end);
        for (i = 0; i < sum; i++) {
            x_list0[i] = malloc(2 * sizeof(double));

            x_list0[i][0] = start + i * dx;
            //printf("[0]: %lf\n", x_list0[i][0]);
            x_list0[i][1] = end - (start + i * dx);
            //printf("[1]: %lf\n", x_list0[i][1]);
        }

        x_list0[sum] = malloc(2 * sizeof(double));
        x_list0[sum][1] = dx;
        x_list0[sum][0] = end - dx;

        return sum + 1;
    }

    sum = flash_calculation_generate_x_number(start,
            end, dx, ncomp);
    sum2 = (int)((end - start) / dx);

    //printf("number: %d\n", sum);

    *x_list = malloc(sum * sizeof(**x_list));
    x_list0 = *x_list;
    for (i = 0; i < sum; i++) {
        x_list0[i] = malloc(ncomp * sizeof(double));
    }

    //printf("ncomp: %d\n", ncomp);
    //printf("end: %e, start: %e, dx: %e\n", end, start, dx);
    //printf("sum2: %d\n", sum2);

    /* add list */
    count = 0;
    for (i = 0; i < sum2; i++) {
        int sum3;
        double x0, **x_list_b;

        x0 = start + i * dx;
        //printf("strat: %lf, x0: %lf, end: %lf, end-x0: %lf\n", 
        //        start, x0, end, end-x0);
        sum3 = flash_calculation_generate_x(start, end - x0, dx, 
                ncomp - 1, min_x, &x_list_b);
        //printf("sum3: %d\n", sum3);

        for (j = 0; j < sum3; j++) {
            x_list0[count][0] = x0;

            for (k = 0; k < ncomp - 1; k++) {
                x_list0[count][k + 1] = x_list_b[j][k];
            }

            count += 1;
            free(x_list_b[j]);
        }

        free(x_list_b);
    }

    while(1) {
        if (end - min_x * (ncomp - 1) > min_x) {
            break;
        }
        else {
            min_x *= 0.5;
        }
    }

    for (i = 0; i < ncomp - 1; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                x_list0[count][j] = end - min_x * (ncomp - 1);
            }
            else {
                x_list0[count][j] = min_x;
            }
        }
        count += 1;
    }

    if (count != sum)
        printf("JAKDFJSLKDJFSLKDJFSDLKFJLSDKJFSDKJFSLDKJFSDLKFJD\n");

    return sum;
}

void flash_calculation_generate_stability_analysis_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, FLASH_STAB_ANN *fsa, char *output)
{
    int i, j, ncomp = comp_list->ncomp;
    STABILITY_MAP *sm;
    double *x;

    x = malloc(ncomp * sizeof(*x));

    printf("   -- Generating stability analysis data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);
        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        sm = flash_calculation_draw_stability_analysis_map(comp_list, x, 
                T_min, T_max, P_min, P_max, dT, dP, fsa, output);

        flash_calculation_stability_map_free(&sm);

        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");
    }
    printf("   -- Done\n");

    free(x);
}

void flash_calculation_generate_stability_analysis_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double P_min, double P_max, double dP, 
        double dxx, FLASH_STAB_ANN *fsa, char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    STABILITY_PM_MAP *sm0, *sm;
    double *x;
    double max_P = 0.0, min_P = 1e10;
    double max_P0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));

    sm = malloc(sizeof(*sm));
    sm->n_unstable = 0;
    sm->unstable_pres = malloc(sizeof(*(sm->unstable_pres)));
    sm->unstable_x = malloc(sizeof(*(sm->unstable_x)));
    sm->n_liquid = 0;
    sm->liquid_pres = malloc(sizeof(*(sm->liquid_pres)));
    sm->liquid_x = malloc(sizeof(*(sm->liquid_x)));
    sm->n_vapor = 0;
    sm->vapor_pres = malloc(sizeof(*(sm->vapor_pres)));
    sm->vapor_x = malloc(sizeof(*(sm->vapor_x)));

    printf("   -- Generating stability analysis data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);
        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        sm0 = flash_calculation_draw_stability_analysis_map_PM(comp_list, x, 
                T, P_min, P_max, dP, 0, dxx, fsa, output);

        sm->unstable_pres = realloc(sm->unstable_pres,
                (sm->n_unstable + sm0->n_unstable) * sizeof(*(sm->unstable_pres)));
        sm->liquid_pres = realloc(sm->liquid_pres,
                (sm->n_liquid + sm0->n_liquid) * sizeof(*(sm->liquid_pres)));
        sm->vapor_pres = realloc(sm->vapor_pres,
                (sm->n_vapor + sm0->n_vapor) * sizeof(*(sm->vapor_pres)));
        sm->unstable_x = realloc(sm->unstable_x,
                (sm->n_unstable + sm0->n_unstable) * sizeof(*(sm->unstable_x)));
        sm->liquid_x = realloc(sm->liquid_x,
                (sm->n_liquid + sm0->n_liquid) * sizeof(*(sm->liquid_x)));
        sm->vapor_x = realloc(sm->vapor_x,
                (sm->n_vapor + sm0->n_vapor) * sizeof(*(sm->vapor_x)));

        for (j = sm->n_unstable; j < sm->n_unstable + sm0->n_unstable; j++) {
            sm->unstable_pres[j] = sm0->unstable_pres[j - sm->n_unstable];

            sm->unstable_x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->unstable_x[j][k] = sm0->unstable_x[j - sm->n_unstable][k];
            }

            if (max_P < sm->unstable_pres[j]) {
                max_P = sm->unstable_pres[j];
            }

            if (min_P > sm->unstable_pres[j]) {
                min_P = sm->unstable_pres[j];
            }
        }
        sm->n_unstable += sm0->n_unstable;

        for (j = sm->n_liquid; j < sm->n_liquid + sm0->n_liquid; j++) {
            sm->liquid_pres[j] = sm0->liquid_pres[j - sm->n_liquid];

            sm->liquid_x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->liquid_x[j][k] = sm0->liquid_x[j - sm->n_liquid][k];
            }

            if (max_P < sm->liquid_pres[j]) {
                max_P = sm->liquid_pres[j];
            }

            if (min_P > sm->liquid_pres[j]) {
                min_P = sm->liquid_pres[j];
            }
        }
        sm->n_liquid += sm0->n_liquid;

        for (j = sm->n_vapor; j < sm->n_vapor + sm0->n_vapor; j++) {
            sm->vapor_pres[j] = sm0->vapor_pres[j - sm->n_vapor];

            sm->vapor_x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->vapor_x[j][k] = sm0->vapor_x[j - sm->n_vapor][k];
            }

            if (max_P < sm->vapor_pres[j]) {
                max_P = sm->vapor_pres[j];
            }

            if (min_P > sm->vapor_pres[j]) {
                min_P = sm->vapor_pres[j];
            }
        }
        sm->n_vapor += sm0->n_vapor;

        flash_calculation_stability_PM_map_free(&sm0);

        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");
    }
    printf("   -- Done\n");

    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
    }

    for (i = 0; i < sm->n_unstable; i++) {
        sm->unstable_pres[i] = (sm->unstable_pres[i] - min_P0) / (max_P0 - min_P0);
    }
    for (i = 0; i < sm->n_liquid; i++) {
        sm->liquid_pres[i] = (sm->liquid_pres[i] - min_P0) / (max_P0 - min_P0);
    }
    for (i = 0; i < sm->n_vapor; i++) {
        sm->vapor_pres[i] = (sm->vapor_pres[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_output_stability_analysis_map_PM(sm, 
            NULL, ncomp, file_scaled);

    free(sm->unstable_pres);
    for (i = 0; i < sm->n_unstable; i++) {
        free(sm->unstable_x[i]);
    }
    free(sm->unstable_x);

    free(sm->liquid_pres);
    for (i = 0; i < sm->n_liquid; i++) {
        free(sm->liquid_x[i]);
    }
    free(sm->liquid_x);

    free(sm->vapor_pres);
    for (i = 0; i < sm->n_vapor; i++) {
        free(sm->vapor_x[i]);
    }
    free(sm->vapor_x);

    free(sm);
    free(x);
}

void flash_calculation_generate_split_calculation_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, FLASH_SPLIT_ANN *fsa, char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    SPLIT_MAP *sm0, *sm;
    double *x;
    double max_T = 0.0, max_P = 0.0, 
           min_T = 1e10, min_P = 1e10;
    double max_T0, max_P0, min_T0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));

    sm = malloc(sizeof(*sm));
    sm->n = 0;
    sm->temp = malloc(sizeof(*(sm->temp)));
    sm->pres = malloc(sizeof(*(sm->pres)));
    sm->F = malloc(sizeof(*(sm->F)));
    sm->K = malloc(sizeof(*(sm->K)));
    sm->x = malloc(sizeof(*(sm->x)));

    printf("   -- Generating split calculation data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);

        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        sm0 = flash_calculation_draw_split_calculation_map(comp_list, x, 
                T_min, T_max, P_min, P_max, dT, dP, fsa, output);

        sm->temp = realloc(sm->temp, 
                (sm->n + sm0->n) * sizeof(*(sm->temp)));
        sm->pres = realloc(sm->pres, 
                (sm->n + sm0->n) * sizeof(*(sm->pres)));
        sm->F = realloc(sm->F, 
                (sm->n + sm0->n) * sizeof(*(sm->F)));
        sm->K = realloc(sm->K, 
                (sm->n + sm0->n) * sizeof(*(sm->K)));
        sm->x = realloc(sm->x, 
                (sm->n + sm0->n) * sizeof(*(sm->x)));

        for (j = sm->n; j < sm->n + sm0->n; j++) {
            sm->temp[j] = sm0->temp[j - sm->n];
            sm->pres[j] = sm0->pres[j - sm->n];
            sm->F[j] = sm0->F[j - sm->n];

            sm->K[j] = malloc(ncomp * sizeof(double));
            sm->x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->K[j][k] = sm0->K[j - sm->n][k];
                sm->x[j][k] = x[k];
            }

            if (max_T < sm->temp[j]) {
                max_T = sm->temp[j];
            }
            if (max_P < sm->pres[j]) {
                max_P = sm->pres[j];
            }

            if (min_T > sm->temp[j]) {
                min_T = sm->temp[j];
            }
            if (min_P > sm->pres[j]) {
                min_P = sm->pres[j];
            }
        }
        sm->n += sm0->n;

        flash_calculation_split_map_free(&sm0);
        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");

    }
    printf("   -- Done\n");

    sprintf(file_scaled, "%s-non-scaled", output);
    flash_calculation_output_split_calculation_map(sm, NULL, 
            ncomp, 1, file_scaled);

    MPI_Allreduce(&max_T, &max_T0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_T, &min_T0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
        printf("Minimal Temperature: %e, Maximal Temperature: %e\n", 
                min_T0, max_T0);
    }

    for (i = 0; i < sm->n; i++) {
        sm->temp[i] = (sm->temp[i] - min_T0) / (max_T0 - min_T0);
        sm->pres[i] = (sm->pres[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_output_split_calculation_map(sm, NULL, 
            ncomp, 1, file_scaled);

    free(sm->pres);
    free(sm->temp);
    free(sm->F);
    for (i = 0; i < sm->n; i++) {
        free(sm->x[i]);
        free(sm->K[i]);
    }
    free(sm->x);
    free(sm->K);
    free(sm);

    free(x);
}

void flash_calculation_generate_split_calculation_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double P_min, double P_max, double dP, 
        double dxx, FLASH_SPLIT_ANN *fsa, char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    SPLIT_PM_MAP *sm0, *sm;
    double *x;
    double max_P = 0.0, min_P = 1e10;
    double max_P0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));

    sm = malloc(sizeof(*sm));
    sm->n = 0;
    sm->pres = malloc(sizeof(*(sm->pres)));
    sm->F = malloc(sizeof(*(sm->F)));
    sm->K = malloc(sizeof(*(sm->K)));
    sm->x = malloc(sizeof(*(sm->x)));

    printf("   -- Generating split calculation data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);

        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        sm0 = flash_calculation_draw_split_calculation_map_PM(comp_list, x, 
                T, P_min, P_max, dP, 0, dxx, fsa, output);

        sm->pres = realloc(sm->pres, 
                (sm->n + sm0->n) * sizeof(*(sm->pres)));
        sm->F = realloc(sm->F, 
                (sm->n + sm0->n) * sizeof(*(sm->F)));
        sm->K = realloc(sm->K, 
                (sm->n + sm0->n) * sizeof(*(sm->K)));
        sm->x = realloc(sm->x, 
                (sm->n + sm0->n) * sizeof(*(sm->x)));

        for (j = sm->n; j < sm->n + sm0->n; j++) {
            sm->pres[j] = sm0->pres[j - sm->n];
            sm->F[j] = sm0->F[j - sm->n];

            sm->K[j] = malloc(ncomp * sizeof(double));
            sm->x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->K[j][k] = sm0->K[j - sm->n][k];
                sm->x[j][k] = sm0->x[j - sm->n][k];
            }

            if (max_P < sm->pres[j]) {
                max_P = sm->pres[j];
            }

            if (min_P > sm->pres[j]) {
                min_P = sm->pres[j];
            }
        }
        sm->n += sm0->n;

        flash_calculation_split_PM_map_free(&sm0);
        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");

    }
    printf("   -- Done\n");

    sprintf(file_scaled, "%s-non-scaled", output);
    flash_calculation_output_split_calculation_map_PM(sm, NULL, 
            ncomp, 1, file_scaled);

    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
    }

    for (i = 0; i < sm->n; i++) {
        sm->pres[i] = (sm->pres[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_output_split_calculation_map_PM(sm, NULL, 
            ncomp, 1, file_scaled);

    free(sm->pres);
    free(sm->F);
    for (i = 0; i < sm->n; i++) {
        free(sm->x[i]);
        free(sm->K[i]);
    }
    free(sm->x);
    free(sm->K);
    free(sm);

    free(x);
}

void flash_calculation_generate_phase_envelope_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    PHASE_DIAGRAM *pd;
    double *x;
    CRITICAL_POINT *cp;
    PHASE_ENVELOPE *pe;
    double max_T = 0.0, max_P = 0.0, 
           min_T = 1e10, min_P = 1e10;
    double max_T0, max_P0, min_T0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));
    pe = malloc(sizeof(*pe));
    pe->n = 0;
    pe->Ps = malloc(sizeof(*(pe->Ps)));
    pe->Ts = malloc(sizeof(*(pe->Ts)));
    pe->xs = malloc(sizeof(*(pe->xs)));

    printf("   -- Generating phase envelope data:\n");

    for (i = 0; i < nx; i++) {
        cp = malloc(sizeof(*cp));
        cp->Tc = 0.0;
        cp->Pc = 0.0;

        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);
        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        pd = flash_calculation_draw_phase_diagram(comp_list, x,
                T_min, T_max, P_min, P_max, NULL, 0, dT, dP, 
                NULL, cp, output);

        pe->Ps = realloc(pe->Ps, (pe->n + pd->pe->n) * sizeof(*(pe->Ps)));
        pe->Ts = realloc(pe->Ts, (pe->n + pd->pe->n) * sizeof(*(pe->Ts)));
        pe->xs = realloc(pe->xs, (pe->n + pd->pe->n) * sizeof(*(pe->xs)));

        for (j = pe->n; j < pe->n + pd->pe->n; j++) {
            pe->Ps[j] = pd->pe->Ps[j - pe->n];
            pe->Ts[j] = pd->pe->Ts[j - pe->n];
            pe->xs[j] = malloc(ncomp * sizeof(*(pe->xs[j])));
            for (k = 0; k < ncomp; k++) {
                pe->xs[j][k] = x[k];
            }

            if (max_T < pe->Ts[j]) {
                max_T = pe->Ts[j];
            }
            if (max_P < pe->Ps[j]) {
                max_P = pe->Ps[j];
            }

            if (min_T > pe->Ts[j]) {
                min_T = pe->Ts[j];
            }
            if (min_P > pe->Ps[j]) {
                min_P = pe->Ps[j];
            }
        }
        pe->n += pd->pe->n;

        flash_calculation_phase_diagram_free(&pd);

        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");
    }
    printf("   -- Done\n");

    MPI_Allreduce(&max_T, &max_T0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_T, &min_T0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
        printf("Minimal Temperature: %e, Maximal Temperature: %e\n", 
                min_T0, max_T0);
    }

    for (i = 0; i < pe->n; i++) {
        pe->Ts[i] = (pe->Ts[i] - min_T0) / (max_T0 - min_T0);
        pe->Ps[i] = (pe->Ps[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_phase_envelope_output(pe, NULL, ncomp, file_scaled);

    free(pe->Ts);
    free(pe->Ps);
    for (i = 0; i < pe->n; i++) {
        free(pe->xs[i]);
    }
    free(pe->xs);
    free(pe);

    free(x);
}

void flash_calculation_generate_phase_envelope_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double dP, double dxx, 
        char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    PHASE_ENVELOPE_PM *pe_pm0, *pe_pm;
    double *x;
    double max_P = 0.0, min_P = 1e10, max_P0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));

    pe_pm = malloc(sizeof(*pe_pm));
    pe_pm->n = 0;
    pe_pm->Ps = malloc(sizeof(*(pe_pm->Ps)));
    pe_pm->xs = malloc(sizeof(*(pe_pm->xs)));

    printf("   -- Generating phase envelope data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);
        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        //for (j = 0; j < ncomp; j++) {
        pe_pm0 = flash_calculation_phase_saturation_envelope_construction_PM(comp_list, 
                x, T, 100.0, dP, 0, dxx, output);

        pe_pm->Ps = realloc(pe_pm->Ps, 
                (pe_pm->n + pe_pm0->n) * sizeof(*(pe_pm->Ps)));
        pe_pm->xs = realloc(pe_pm->xs, 
                (pe_pm->n + pe_pm0->n) * sizeof(*(pe_pm->xs)));

        for (j = pe_pm->n; j < pe_pm->n + pe_pm0->n; j++) {
            pe_pm->Ps[j] = pe_pm0->Ps[j - pe_pm->n];
            pe_pm->xs[j] = malloc(ncomp * sizeof(*(pe_pm->xs[j])));
            for (k = 0; k < ncomp; k++) {
                pe_pm->xs[j][k] = pe_pm0->xs[j - pe_pm->n][k];
            }

            if (max_P < pe_pm->Ps[j]) {
                max_P = pe_pm->Ps[j];
            }
            if (min_P > pe_pm->Ps[j]) {
                min_P = pe_pm->Ps[j];
            }
        }
        pe_pm->n += pe_pm0->n;

        flash_calculation_phase_envelope_pm_free(&pe_pm0);
        //}

        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");
    }
    printf("   -- Done\n");

    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
    }

    for (i = 0; i < pe_pm->n; i++) {
        pe_pm->Ps[i] = (pe_pm->Ps[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_phase_envelope_PM_output(pe_pm, ncomp, 0, file_scaled);

    free(pe_pm->Ps);
    for (i = 0; i < pe_pm->n; i++) {
        free(pe_pm->xs[i]);
    }
    free(pe_pm->xs);
    free(pe_pm);

    free(x);
}


