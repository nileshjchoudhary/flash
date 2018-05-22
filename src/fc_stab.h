#ifndef FC_STAB
typedef struct STABILITY_MAP_ {
    int n_unstable;
    double *unstable_pres;
    double *unstable_temp;

    int n_liquid;
    double *liquid_pres;
    double *liquid_temp;

    int n_vapor;
    double *vapor_pres;
    double *vapor_temp;

} STABILITY_MAP;

typedef struct STABILITY_PM_MAP_ {
    int n_unstable;
    double *unstable_pres;
    double **unstable_x;

    int n_liquid;
    double *liquid_pres;
    double **liquid_x;

    int n_vapor;
    double *vapor_pres;
    double **vapor_x;

} STABILITY_PM_MAP;

double * flash_calculation_estimate_K(EOS *eos, double *K);
double * flash_calculation_stability_analysis_initial_estimate(PHASE *phase);
void flash_calculation_calculate_trial_phase_composition(double *X_t, double *x, int ncomp);
void flash_calculation_SS_method_update_X(PHASE *phase, PHASE *phase_t, double *X_t);
void flash_calculation_calculate_trial_phase_composition_derivative(double *X_t, double *dx_t, int ncomp);
void flash_calculation_calculate_stability_equilibrium_equation(PHASE *phase, PHASE *phase_t, 
        double *X_t, double *D);
void flash_calculation_calculate_stability_equilibrium_equation_derivative(PHASE *phase_t, double *dx_t, 
        double *X_t, double *dD);
double flash_calculation_calculate_stability_residual(PHASE *phase, PHASE *phase_t, double *X_t, double *res);
void flash_calculation_QNSS_method_update_X(double *dD, double *D, double *X_t, int ncomp);
int flash_calculation_check_stability(double *X_t, double *z, int ncomp)
int flash_calculation_stability_analysis_QNSS(PHASE *phase, double *K, double tol)


void flash_calculation_output_stability_analysis_map(STABILITY_MAP *sm, 
        double *comp_X, int ncomp, char *output_name);
void flash_calculation_output_stability_analysis_map_PM(STABILITY_PM_MAP *sm, 
        double *comp_X, int ncomp, char *output_name);
STABILITY_MAP * flash_calculation_draw_stability_analysis_map(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double dT, double dP, FLASH_STAB_ANN *fsa, char *output_name);
STABILITY_PM_MAP * flash_calculation_draw_stability_analysis_map_PM(COMP_LIST *comp_list, 
        double *comp_X, double T, double P_min, double P_max, double dP, int selected_component, 
        double dxx, FLASH_STAB_ANN *fsa, char *output_name);
void flash_calculation_stability_map_free(STABILITY_MAP **sm);
void flash_calculation_stability_PM_map_free(STABILITY_PM_MAP **sm);
#define FC_STAB
#endif
