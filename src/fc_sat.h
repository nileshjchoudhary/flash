#ifndef FC_SAT
void flash_calculation_saturation_calculation_search_Pu_Ps(EOS *eos, double *z, double T, double P_est, 
	int search, double search_step, int Ps_found, int Pu_found, double *Pu_Ps, double P_max);
void flash_calculation_saturation_calculation_search_Tu_Ts(EOS *eos, double *z, double P, double T_est, int search, 
	double search_step, int Ts_found, int Tu_found, double *Tu_Ts);
double flash_calculation_saturation_calculation_search_boundary_pressure(EOS *eos, double *z,
	double Pu, double Ps, int search, double tol, int search_status, double *K);
void flash_calculation_saturation_calculation_calculate_X(double *K, double *z, double *X, int ncomp);
void flash_calculation_aturation_calculation_calculate_X_derivative(double *z, double *dX_dK, int ncomp);
void flash_calculation_saturation_calculation_calculate_fugacity_ratio(PHASE *phase, PHASE *phase_t, 
	double *Y, double *Ri);
void flash_calculation_saturation_calculation_update_X_with_fugacity_ratio(double *Y, double *Ri, int ncomp);
double flash_calculation_saturation_calculation_calculate_Q_value(double *Y, int ncomp);
double flash_calculation_saturation_calculation_calculate_Q_derivative(PHASE *phase, PHASE *phase_t, 
	double *Y, double *Ri);
double flash_calculation_saturation_calculation(EOS *eos, double *z, double T, double P_est, 
        int search, double search_step, double P_max);
#define FC_SAT
#endif
