#include "fc.h"

COMP_LIST * flash_calculation_component_new(char *prop_file, char *bin_file)
{
    FILE *fn;
    char line[10240];
    CSV_LINE **csv_value, *csv_line0; 
    int nline, i, j, ncomp;
    COMP_LIST *comp_list;
    COMP *comp;

    csv_value = malloc(20 * sizeof(*csv_value));

    fn = fopen(prop_file, "r");
    nline = 0;
    while (fgets(line, 10240, fn)){
        csv_line0 = get_CSV_line(line, NULL);
        csv_value[nline] = csv_line0;

        nline++;
    }
    fclose(fn);

    ncomp = csv_value[0]->n - 1;
    comp = malloc(ncomp * sizeof(*comp));

    for (i = 0; i < nline; i++) {
        char *name;

        name = csv_value[i]->value[0];

        if (strcmp(name, "PC") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].PC = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "TC") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].TC = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "VC") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].VC = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "AC") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].AC = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "MW") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].MW = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "VS") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].VS = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "HYD") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].HYD = atoi(csv_value[i]->value[j + 1]);
            }
        }
    }

    for (i = 0; i < nline; i++) {
        free_CSV_line(&(csv_value[i]));
    }

    fn = fopen(bin_file, "r");
    nline = 0;
    while (fgets(line, 10240, fn)){
        csv_line0 = get_CSV_line(line, NULL);
        csv_value[nline] = csv_line0;

        nline++;
    }
    fclose(fn);

    for (i = 0; i < nline; i++) {
        comp[i].binary = malloc(ncomp * sizeof(*(comp[i].binary)));

        for (j = 0; j < ncomp; j++) {    
            comp[i].binary[j] = atof(csv_value[i]->value[j]);
        }
    }

    for (i = 0; i < nline; i++) {
        free_CSV_line(&(csv_value[i]));
    }

    comp_list = malloc(sizeof(*comp_list));
    comp_list->ncomp = ncomp;
    comp_list->comp = comp;

    free(csv_value);

    return comp_list;
}

void flash_calculation_component_free(COMP_LIST **comp_list)
{
    int i;
    COMP_LIST *comp_list0 = *comp_list;

    for (i = 0; i < comp_list0->ncomp; i++) {
        free(comp_list0->comp[i].binary);
    }

    free(comp_list0->comp);

    free(*comp_list);
}

double * flash_calculation_composition_new(char *z_file)
{
    FILE *fn;
    char line[10240];
    CSV_LINE **csv_value, *csv_line0; 
    int nline, i;
    double *z;

    csv_value = malloc(20 * sizeof(*csv_value));

    fn = fopen(z_file, "r");
    nline = 0;
    while (fgets(line, 10240, fn)){
        csv_line0 = get_CSV_line(line, NULL);
        csv_value[nline] = csv_line0;

        nline++;
    }
    fclose(fn);

    z = malloc(nline * sizeof(*z));

    for (i = 0; i < nline; i++) {
        z[i] = atof(csv_value[i]->value[0]);
    }

    for (i = 0; i < nline; i++) {
        free_CSV_line(&(csv_value[i]));
    }

    free(csv_value);

    return z;
}
