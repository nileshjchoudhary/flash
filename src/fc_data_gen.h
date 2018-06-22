#ifndef FC_DATA_GEN
int flash_calculation_generate_x(double start, double end, double dx, 
        int ncomp, double min_x, double ***x_list);
int flash_calculation_generate_x_new(int mole, int *mole_range, int *mole_dx,
        int ncomp, double ***x_list);
int flash_calculation_generate_x_new_2(int mole, int *mole_range, int *mole_d,
        int ncomp, double ***x_list);
void flash_calculation_generate_stability_analysis_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, FLASH_STAB_ANN *fsa, char *output);
void flash_calculation_generate_stability_analysis_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double P_min, double P_max, double dP, 
        double dxx, FLASH_STAB_ANN *fsa, char *output);
void flash_calculation_generate_split_calculation_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, FLASH_SPLIT_ANN *fsa, char *output);
void flash_calculation_generate_split_calculation_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double P_min, double P_max, double dP, 
        double *comp_range, double dxx, FLASH_SPLIT_ANN *fsa, char *output);
void flash_calculation_generate_phase_envelope_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, char *output);
void flash_calculation_generate_phase_envelope_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double dP, double *comp_range, double dxx, 
        double P_max, char *output);
#define FC_DATA_GEN
#endif
