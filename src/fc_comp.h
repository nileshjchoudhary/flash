#ifndef FC_COMP
typedef struct COMP_ {
    char name[20];
    double PC;
    double TC;
    double VC;
    double AC;
    double MW;
    double VS;
    int HYD;
    double *binary;
} COMP;

typedef struct COMP_LIST_ {
    int ncomp;
    COMP *comp;

} COMP_LIST;

COMP_LIST * flash_calculation_component_new(char *prop_file, char *bin_file);
void flash_calculation_component_free(COMP_LIST **comp_list);
double * flash_calculation_composition_new(char *z_file);

#define FC_COMP
#endif
