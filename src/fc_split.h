#ifndef FC_SPLIT
typedef struct SPLIT_MAP_ {
    int n;
    double *temp;
    double *pres;
    double *F;
    double **K;
    double **x;

} SPLIT_MAP;

typedef struct SPLIT_PM_MAP_ {
    int n;
    double *pres;
    double *F;
    double **K;
    double **x;

} SPLIT_PM_MAP;

void flash_calculation_calculate_composition(double *K, double *z, double F_v, 
        double *x_l, double *x_v, int ncomp);
void flash_calculation_calculate_composition_derivative(double *K, double *z, double F_v, 
        double *dx_l, double *dx_v, int ncomp);
double flash_calculation_calculate_equilibrium_equation(PHASE *phase_L, PHASE *phase_V, double *G);
void flash_calculation_calculate_equilibrium_equation_derivative(PHASE *phase_L, PHASE *phase_V, 
        double *dx_l, double *dx_v, double *dG);
void flash_calculation_QNSS_method_update_K(double *dG, double *G, double *K, int ncomp);
double flash_calculation_calculate_RachfordRice_equation_value(double *K, double *z, double n_V, int ncomp);
double flash_calculation_calculate_RachfordRice_equation_derivative(double *K, double *z, double n_V, int ncomp);
double flash_calculation_solve_RachfordRice_equation(double *K, double *z, double n_V0, int ncomp);
void flash_calculation_SS_method_update_K(double *fug_L, double *fug_V, double *K, int ncomp);
double flash_calculation_two_phase_flash_calculation_calculate_initial_F(double *K, double *z, int ncomp);
double flash_calculation_two_phase_flash_Calculation_QNSS(EOS *eos, double *z, 
        double *K, double Fv, double tol);


void flash_calculation_output_split_calculation_map(SPLIT_MAP *sm, 
        double *comp_X, int ncomp, int filter, char *output_name);
void flash_calculation_output_split_calculation_map_PM(SPLIT_PM_MAP *sm, 
        double *comp_X, int ncomp, int filter, char *output_name);
SPLIT_MAP * flash_calculation_draw_split_calculation_map(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double dT, double dP, FLASH_SPLIT_ANN *fsa, char *output_name);
SPLIT_PM_MAP * flash_calculation_draw_split_calculation_map_PM(COMP_LIST *comp_list, 
        double *comp_X, double T, double P_min, double P_max, double dP, 
        int selected_component, double dxx, FLASH_SPLIT_ANN *fsa, char *output_name);
void flash_calculation_split_map_free(SPLIT_MAP **sm);
void flash_calculation_split_PM_map_free(SPLIT_PM_MAP **sm);

double flash_calculation_split_time_cost(void);
int flash_calculation_split_iteration_number(void);
int flash_calculation_split_failure_number(void);
double flash_calculation_split_pred_time_cost(void);

#define FC_SPLIT
#endif
