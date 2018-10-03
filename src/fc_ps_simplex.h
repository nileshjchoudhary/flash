#ifndef FC_PS_SIMPLEX

typedef struct PS_SIMPLEX_ISOTHERM_ {
    double T;

    int n;
    double *Ps_u;
    double *Ps_l;
    double **xs;

} PS_SIMPLEX_ISOTHERM;

typedef struct SPLIT_SIMPLEX_ISOTHERM_ {
    double T;
    
    int n;
    int *nP;
    double **P;
    double **Fv;
    double ***K;
    double **xs;

} SPLIT_SIMPLEX_ISOTHERM;

typedef struct PS_SIMPLEX_ {
    int nT;
    PS_SIMPLEX_ISOTHERM *ps_isotherm;

} PS_SIMPLEX;

typedef struct SET_NO_LIST_ {
    int set_begin;
    int set_end;
    int set_init;

    struct SET_NO_LIST_ *parent;
    struct SET_NO_LIST_ *next;
    struct SET_NO_LIST_ *previous;

    int nchild;
    struct SET_NO_LIST_ **child;

} SET_NO_LIST;

int flash_calculation_generate_simplex(double start, double end, 
        double dx, int ncomp, double ***x_list);
PS_SIMPLEX_ISOTHERM * 
flash_calculation_saturation_pressure_simplex_isotherm(COMP_LIST *comp_list,
        double **z, int nz, double *z_range, SET_NO_LIST *set_no_list, double T, 
        double Ps_u_est, double Ps_l_est, double dP, double P_max);
void
flash_calculation_saturation_pressure_simplex_isotherm_free(PS_SIMPLEX_ISOTHERM **ps);
void flash_calculation_saturation_pressure_simplex_isotherm_output(PS_SIMPLEX_ISOTHERM *ps, 
        int ncomp, char *output);
SET_NO_LIST * flash_calculation_generate_simplex_set_no(double start, 
        double end, double dx, int ncomp);
void flash_calculation_generate_simplex_set_no_free(SET_NO_LIST **set_no_list);
void flash_calculation_generate_simplex_set_no_print(SET_NO_LIST *set_no_list);
PS_SIMPLEX_ISOTHERM *
flash_calculation_saturation_pressure_simplex_isotherm_data(COMP_LIST *comp_list,
    double T, double dx, double *z_range, double dP, double P_max, char *output);

SPLIT_SIMPLEX_ISOTHERM *
flash_calculation_split_simplex_isotherm(COMP_LIST *comp_list,
        PS_SIMPLEX_ISOTHERM *ps, double dP);
void flash_calculation_split_simplex_isotherm_output(SPLIT_SIMPLEX_ISOTHERM *sp, 
        int ncomp, char *output);
void flash_calculation_split_simplex_isotherm_free(SPLIT_SIMPLEX_ISOTHERM **sp);
SPLIT_SIMPLEX_ISOTHERM *
flash_calculation_split_simplex_isotherm_data(COMP_LIST *comp_list,
    PS_SIMPLEX_ISOTHERM *ps, double dP, char *output);
void flash_calculation_simplex_isotherm_data(COMP_LIST *comp_list, 
        double T, double dx, double *z_range, double dP, double P_max, char *output);

#define FC_PS_SIMPLEX
#endif
