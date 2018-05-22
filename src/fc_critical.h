#ifndef FC_CRITICAL
typedef struct CRITICAL_POINT_ {
    double Pc;
    double Tc;
} CRITICAL_POINT;

double flash_calculation_calculate_initial_K_derivative_with_pressure(EOS *eos);
double flash_calculation_calculate_initial_K_derivative_with_temperature(EOS *eos);
double * flash_calculation_calculate_status_at_fixed_temperature_and_Fv(EOS *eos, double *z, 
        double T, double Fv, double P_max, double P_min, double P_est, double *K);
double * flash_calculation_calculate_status_at_fixed_pressure_and_Fv(EOS *eos, double *z, 
        double P, double Fv, double T_max, double T_min, double T_est, double *K); 
double flash_calculation_critial_point_calculate_init_T_guess(EOS *eos, double *z);
void flash_calculation_critical_point_volume_functions(double kappa, double sigma1, double sigma2, 
	double *F);
void flash_calculation_critical_point_calculate_alpha_beta(PHASE *phase, double *alpha, double *beta);
void flash_calculation_critical_point_calculate_Q_matrix(double *F, PHASE *phase, double *alpha, double *beta, double *Q);
double flash_calculation_calculate_matrix_determinant(double *Q, int dim);
double flash_calculation_critical_point_calculate_Q_det(double *F, PHASE *phase, 
        double Tc, double *Qij);
double flash_calculation_critical_point_calculate_Q_det_derivative(double *F, 
        PHASE *phase, double Tc);
double flash_calculation_critical_point_calculate_cubic_form(PHASE *phase, 
        double *N, double *alpha, double *beta, double kappa);
double flash_calculation_critical_point_calculate_cubic_form_derivative(PHASE *phase, 
        double *N, double *alpha, double *beta, double kappa);
void flash_calculation_critical_point_calculate_normal_N(double *Q, int ncomp, 
        double *x);
CRITICAL_POINT * flash_calculation_critical_point_calculation(EOS *eos, 
        double *z, double kappa_init, double Tc_init);
void flash_calculation_critical_point_free(CRITICAL_POINT **cp);
#define FC_CRITICAL
#endif
