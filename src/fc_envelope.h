#ifndef PC_ENVELOPE
typedef struct PHASE_ENVELOPE_ {
    int n;

    double *Ps;
    double *Ts;
    double **xs;

} PHASE_ENVELOPE;

typedef struct PHASE_ENVELOPE_PM_ {
    int n;
    double *Ps;
    double **xs;

} PHASE_ENVELOPE_PM;
double flash_calculation_phase_diagram_construction_calculate_equations(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *G);
void flash_calculation_phase_diagram_construction_calculate_equations_derivative(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *dG);
void flash_calculation_phase_diagram_construction_calculate_equations_derivative_with_K_T(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *dG);
double * flash_calculation_phase_diagram_construction_update_variables(double *dG, double *G, int dim);
int flash_calculation_search_unstable_temperature(COMP_LIST *comp_list, 
        double *comp_X, double P, double *T_list, int nT);
int flash_calculation_search_stable_temperature(COMP_LIST *comp_list, 
        double *comp_X, double P, double *T_list, int nT);
PHASE_ENVELOPE * flash_calculation_phase_saturation_envelope_construction(EOS *eos, 
        double *z, double T_start, double T_end, double dT, double P_est, double dP, double P_max);
void flash_calculation_phase_envelope_PM_output(PHASE_ENVELOPE_PM *pe_pm,
        int ncomp, int selected_component, char *output);
PHASE_ENVELOPE_PM * flash_calculation_phase_saturation_envelope_construction_PM(COMP_LIST *comp_list, 
        double *z, double T, double P_est, double dP, int selected_component, double dx, 
        char *output);

void flash_calculation_phase_envelope_output(PHASE_ENVELOPE *pe,
        double *comp_X, int ncomp, char *output_name);
void flash_calculation_phase_envelope_free(PHASE_ENVELOPE **pe);
void flash_calculation_phase_envelope_pm_free(PHASE_ENVELOPE_PM **pe_pm);
#define PC_ENVELOPE
#endif
