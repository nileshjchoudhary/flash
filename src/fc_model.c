#include "fc.h"

FLASH_MODEL * flash_calculation_model_new(char *model_file)
{
    FILE *fn;
    char line[10240];
    CSV_LINE **csv_value, *csv_line0; 
    int nline, i, j;
    FLASH_MODEL *fm;

    fm = malloc(sizeof(*fm));
    fm->type = malloc(100 * sizeof(*(fm->type)));
    fm->prop_file = malloc(100 * sizeof(*(fm->prop_file)));
    fm->bin_file = malloc(100 * sizeof(*(fm->bin_file)));
    fm->z_file = malloc(100 * sizeof(*(fm->z_file)));
    fm->output = malloc(100 * sizeof(*(fm->output)));
    fm->nF = 0;
    fm->F = NULL;
    fm->mole_range = NULL;
    fm->mole_dx = NULL;
    fm->mole_d = NULL;
    fm->split_ann = NULL;
    fm->split_ann_level = 0;
    fm->stab_ann = NULL;
    fm->stab_ann_level = 0;
    fm->ann_trans = malloc(100 * sizeof(*(fm->ann_trans)));

    sprintf(fm->ann_trans, "none");

    csv_value = malloc(20 * sizeof(*csv_value));

    fn = fopen(model_file, "r");
    nline = 0;
    while (fgets(line, 10240, fn)){
        csv_line0 = get_CSV_line(line, NULL);
        csv_value[nline] = csv_line0;

        nline++;
    }
    fclose(fn);

    for (i = 0; i < nline; i++) {
        char *name;

        name = csv_value[i]->value[0];

        if (strcmp(name, "type") == 0) {
            strcpy(fm->type, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "prop_file") == 0) {
            strcpy(fm->prop_file, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "bin_file") == 0) {
            strcpy(fm->bin_file, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "z_file") == 0) {
            strcpy(fm->z_file, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "T_min") == 0) {
            fm->T_min = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "T_max") == 0) {
            fm->T_max = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "P_min") == 0) {
            fm->P_min = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "P_max") == 0) {
            fm->P_max = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "dT") == 0) {
            fm->dT = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "dP") == 0) {
            fm->dP = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "mole") == 0) {
            fm->mole = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "dx") == 0) {
            fm->dx = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "F") == 0) {
            fm->nF = csv_value[i]->n - 1;
            fm->F = malloc(fm->nF * sizeof(*(fm->F)));

            for (j = 0; j < fm->nF; j++) {
                fm->F[j] = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "mole_range") == 0) {
            fm->n_mole_range = csv_value[i]->n - 1;
            fm->mole_range = malloc(fm->n_mole_range * sizeof(*(fm->mole_range)));

            for (j = 0; j < fm->n_mole_range; j++) {
                fm->mole_range[j] = atoi(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "mole_dx") == 0) {
            fm->n_mole_dx = csv_value[i]->n - 1;
            fm->mole_dx = malloc(fm->n_mole_dx * sizeof(*(fm->mole_dx)));

            for (j = 0; j < fm->n_mole_dx; j++) {
                fm->mole_dx[j] = atoi(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "mole_d") == 0) {
            fm->n_mole_d = csv_value[i]->n - 1;
            fm->mole_d = malloc(fm->n_mole_d * sizeof(*(fm->mole_d)));

            for (j = 0; j < fm->n_mole_d; j++) {
                fm->mole_d[j] = atoi(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "T") == 0) {
            fm->T = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "selected_component") == 0) {
            fm->selected_component = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "dxx") == 0) {
            fm->dxx = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "ann_trans") == 0) {
            strcpy(fm->ann_trans, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "split_ann") == 0) {
            fm->split_ann = malloc(100 * sizeof(*(fm->split_ann)));

            strcpy(fm->split_ann, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "split_ann_level") == 0) {
            fm->split_ann_level = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "stab_ann") == 0) {
            fm->stab_ann = malloc(100 * sizeof(*(fm->stab_ann)));

            strcpy(fm->stab_ann, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "stab_ann_level") == 0) {
            fm->stab_ann_level = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "stab_ann_safeguard") == 0) {
            fm->stab_ann_safeguard = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "stab_ann_delta_p") == 0) {
            fm->stab_ann_delta_p = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "output") == 0) {
            strcpy(fm->output, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "verb") == 0) {
            fm->verb = atoi(csv_value[i]->value[1]);
        }
    }

    return fm;
}

void flash_calculation_flash_model_free(FLASH_MODEL **fm)
{
    FLASH_MODEL *fm0 = *fm;

    free(fm0->type);
    free(fm0->prop_file);
    free(fm0->bin_file);
    free(fm0->z_file);
    if (fm0->F != NULL) {
        free(fm0->F);
    }
    free(fm0->output);

    if (fm0->split_ann != NULL)
        free(fm0->split_ann);

    free(*fm);
}
