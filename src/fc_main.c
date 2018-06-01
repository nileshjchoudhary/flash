#include "fc.h"

int main(int argc, char **argv) 
{
    char *model_file, *prop_file, *binary_file, *z_file;
    COMP_LIST *comp_list;
    double *x;
    int i, j;
    int nprocs, myrank;
    FLASH_MODEL *fm;
    FLASH_SPLIT_ANN *fsa = NULL;
    FLASH_STAB_ANN *fsta = NULL;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    printf("RANK[%d]: argc: %d, argv[1]: %s\n", myrank, argc, argv[1]);

    model_file = argv[1];
    fm = flash_calculation_model_new(model_file);

    prop_file = fm->prop_file;
    binary_file = fm->bin_file;
    z_file = fm->z_file;

    printf("prop_name: %s\n", prop_file);
    printf("binary_name: %s\n", binary_file);
    printf("z_name: %s\n", z_file);

    comp_list = flash_calculation_component_new(prop_file, binary_file);
    x = flash_calculation_composition_new(z_file);

    if (fm->split_ann != NULL && fm->split_ann_level >= 0) {
        fsa = flash_calculation_split_ann_model_new(fm->split_ann, 
                fm->split_ann_level, comp_list->ncomp);
    }

    if (fm->stab_ann != NULL && fm->stab_ann_level >= 0) {
        fsta = flash_calculation_stab_ann_model_new(fm->stab_ann, 
                fm->stab_ann_level, fm->stab_ann_safeguard, 
                fm->stab_ann_delta_p);
    }

    if (strcmp(fm->type, "stability") == 0) {
        STABILITY_MAP *sm;

        sm = flash_calculation_draw_stability_analysis_map(comp_list, x, 
                fm->T_min, fm->T_max, fm->P_min, fm->P_max, fm->dT, fm->dP, 
                fsta, fm->output);

        flash_calculation_stability_map_free(&sm);
    }

    if (strcmp(fm->type, "stability_PM") == 0) {
        STABILITY_PM_MAP *sm;

        sm = flash_calculation_draw_stability_analysis_map_PM(comp_list, x, 
                fm->T, fm->P_min, fm->P_max, fm->dP, fm->selected_component, 
                fm->dxx, fsta, fm->output);

        flash_calculation_stability_PM_map_free(&sm);
    }

    if (strcmp(fm->type, "stability_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

#if 0
        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);
#endif

        nx = flash_calculation_generate_x_new(fm->mole, fm->mole_range, fm->mole_dx,
                comp_list->ncomp, &x_list);


        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%04d", fm->output, myrank);
        }

        flash_calculation_generate_stability_analysis_data(comp_list, nx_rank, 
                x_list + nx_begin, fm->T_min, fm->T_max, fm->P_min, fm->P_max, 
                fm->dT, fm->dP, fsta, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }

    if (strcmp(fm->type, "stability_PM_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);

        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%04d", fm->output, myrank);
        }

        flash_calculation_generate_stability_analysis_PM_data(comp_list, nx_rank, 
                x_list + nx_begin, fm->T, fm->P_min, fm->P_max, fm->dP, 
                fm->dxx, fsta, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }

    if (strcmp(fm->type, "split") == 0) {
        SPLIT_MAP *sm;
        char file_scaled[100];

        sm = flash_calculation_draw_split_calculation_map(comp_list, x, 
                fm->T_min, fm->T_max, fm->P_min, fm->P_max, fm->dT, fm->dP, 
                fsa, fm->output);

        for (i = 0; i < sm->n; i++) {
            sm->temp[i] = (sm->temp[i] - fm->T_min) / (fm->T_max - fm->dT - fm->T_min);
            sm->pres[i] = (sm->pres[i] - fm->P_min) / (fm->P_max - fm->dP - fm->P_min);
        }

        sprintf(file_scaled, "%s-scaled", fm->output);
        flash_calculation_output_split_calculation_map(sm, x, 
                comp_list->ncomp, 1, file_scaled);

        flash_calculation_split_map_free(&sm);
    }

    if (strcmp(fm->type, "split_PM") == 0) {
        SPLIT_PM_MAP *sm;
        char file_scaled[100];

        sm = flash_calculation_draw_split_calculation_map_PM(comp_list, x, 
                fm->T, fm->P_min, fm->P_max, fm->dP, fm->selected_component,
                fm->dxx, fsa, fm->output);

        for (i = 0; i < sm->n; i++) {
            sm->pres[i] = (sm->pres[i] - fm->P_min) / (fm->P_max - fm->dP - fm->P_min);
        }

        sprintf(file_scaled, "%s-scaled", fm->output);
        flash_calculation_output_split_calculation_map_PM(sm, NULL, 
                comp_list->ncomp, 1, file_scaled);

        flash_calculation_split_PM_map_free(&sm);
    }

    if (strcmp(fm->type, "split_data") == 0) {
        int nx, nx_rank;
        double **x_list, **x_list_rank;
        char output_rank[100];

        nx = flash_calculation_generate_x_new_2(fm->mole, 
                fm->mole_range, fm->mole_d,
                comp_list->ncomp, &x_list);

        nx_rank = nx / nprocs;
        if (myrank < nx - nx_rank * nprocs) {
            nx_rank++;
        }

        x_list_rank = malloc(nx_rank * sizeof(*x_list_rank));
        j = 0;
        for (i = 0; i < nx; i++) {
            if (i % nprocs == myrank) {
                x_list_rank[j] = x_list[i];
                j++;
            }
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%03d", fm->output, myrank);
        }

        flash_calculation_generate_split_calculation_data(comp_list, 
                nx_rank, x_list_rank, fm->T_min, fm->T_max, fm->P_min, 
                fm->P_max, fm->dT, fm->dP, fsa, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }

    if (strcmp(fm->type, "split_PM_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);

        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%04d", fm->output, myrank);
        }

        flash_calculation_generate_split_calculation_PM_data(comp_list, 
                nx_rank, x_list + nx_begin, fm->T, fm->P_min, fm->P_max, 
                fm->dP, fm->dxx, fsa, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }



    if (strcmp(fm->type, "diagram") == 0) {
        PHASE_DIAGRAM *pd;

        pd = flash_calculation_draw_phase_diagram(comp_list, x, fm->T_min, fm->T_max,
                fm->P_min, fm->P_max, fm->F, fm->nF, fm->dT, fm->dP, NULL, NULL,
                fm->output);

        flash_calculation_phase_diagram_free(&pd);
    }

    if (strcmp(fm->type, "envelope_PM") == 0) {
        PHASE_ENVELOPE_PM *pe_pm;

        pe_pm = flash_calculation_phase_saturation_envelope_construction_PM(comp_list, 
                x, fm->T, 100.0, fm->dP, fm->selected_component, NULL, fm->dxx, fm->P_max, 
                fm->output);

        flash_calculation_phase_envelope_pm_free(&pe_pm);
    }

    if (strcmp(fm->type, "envelope_data") == 0) {
        int nx, nx_rank;
        double **x_list, **x_list_rank;
        char output_rank[100];

        nx = flash_calculation_generate_x_new_2(fm->mole, 
                fm->mole_range, fm->mole_d,
                comp_list->ncomp, &x_list);

        nx_rank = nx / nprocs;
        if (myrank < nx - nx_rank * nprocs) {
            nx_rank++;
        }

        x_list_rank = malloc(nx_rank * sizeof(*x_list_rank));
        j = 0;
        for (i = 0; i < nx; i++) {
            if (i % nprocs == myrank) {
                x_list_rank[j] = x_list[i];
                j++;
            }
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%03d", fm->output, myrank);
        }

        flash_calculation_generate_phase_envelope_data(comp_list, 
                nx_rank, x_list_rank, fm->T_min, fm->T_max, fm->P_min, 
                fm->P_max, fm->dT, fm->dP, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
        free(x_list_rank);
    }

    if (strcmp(fm->type, "envelope_PM_data") == 0) {
        int nx, nx_rank;
        double **x_list, **x_list_rank, comp_range[2];
        char output_rank[100];

        nx = flash_calculation_generate_x_new_2(fm->mole, 
                fm->mole_range, fm->mole_d,
                comp_list->ncomp, &x_list);

        if (fm->mole_range == NULL) {
            comp_range[0] = 1.0 / fm->mole;
            comp_range[1] = 1.0;
        }
        else {
            comp_range[0] = (double)fm->mole_range[0] / fm->mole;
            comp_range[1] = (double)fm->mole_range[1] / fm->mole;
        }

        nx_rank = nx / nprocs;
        if (myrank < nx - nx_rank * nprocs) {
            nx_rank++;
        }

        x_list_rank = malloc(nx_rank * sizeof(*x_list_rank));
        j = 0;
        for (i = 0; i < nx; i++) {
            if (i % nprocs == myrank) {
                x_list_rank[j] = x_list[i];
                j++;
            }
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%03d", fm->output, myrank);
        }

        flash_calculation_generate_phase_envelope_PM_data(comp_list, 
                nx_rank, x_list_rank, fm->T, fm->dP, 
                comp_range, fm->dxx, fm->P_max, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
        free(x_list_rank);
    }







#if 0
    if (strcmp(fm->type, "thermal_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list, **x_list_rank;
        char output_rank[100];

        nx = flash_calculation_generate_x_new_2(fm->mole, 
                fm->mole_range, fm->mole_d,
                comp_list->ncomp, &x_list);

        nx_rank = nx / nprocs;
        if (myrank < nx - nx_rank * nprocs) {
            nx_rank++;
        }

        x_list_rank = malloc(nx_rank * sizeof(*x_list_rank));
        j = 0;
        for (i = 0; i < nx; i++) {
            if (i % nprocs == myrank) {
                x_list_rank[j] = x_list[i];
                j++;
            }
        }

        printf("rank[%d]---  nx: %d, nx_rank: %d, nx_begin: %d\n", myrank, nx, nx_rank, nx_begin);

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%03d", fm->output, myrank);
        }

        flash_calculation_generate_thermal_data(comp_list, 
                nx_rank, x_list_rank, fm->T_min, fm->T_max, fm->P_min, 
                fm->P_max, fm->dT, fm->dP, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
        free(x_list_rank);
    }
#endif











    printf("Stability calculation iteration: %d\n", 
            flash_calculation_stability_iteration_number());
    printf("Stability calculation ANN prediction time: %lf\n", 
            flash_calculation_stability_pre_time_cost());
    printf("Stability calculation solve time: %lf\n", 
            flash_calculation_stability_time_cost());

    printf("Split calculation failure: %d\n",
            flash_calculation_split_failure_number());
    printf("Split calculation iteration: %d\n", 
            flash_calculation_split_iteration_number());
    printf("Split calculation ANN prediction time: %lf\n", 
            flash_calculation_split_pred_time_cost());
    printf("Split calculation solve time: %lf\n", 
            flash_calculation_split_time_cost());


    if (fsa != NULL) {
        flash_calculation_split_ann_model_free(&fsa);
    }

    free(x);
    flash_calculation_flash_model_free(&fm);
    flash_calculation_component_free(&comp_list);

    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}
