#include "fc.h"

/* ## 5. Saturation Calculations
# The following code is designed for saturation calculations.
*/

/* ### Search pressure at a given temperature
# At a given temperature, search a pressure $P_s$ at which a single phase exists and 
a pressure $P_u$ at which two phases co-exist. The stability test is used in this function.
*/
void flash_calculation_saturation_calculation_search_Pu_Ps(EOS *eos, double *z, double T, double P_est, 
	int search, double search_step, int Ps_found, int Pu_found, double *Pu_Ps, double P_max)
{
	/* search: 0 upper
			   1 down */

	int ncomp = eos->ncomp, status;
	double *K, Ps, Pu;
    PHASE *phase;

	eos->temp = T;
    eos->pres = P_est;
    phase = flash_calculation_phase_new(eos, z);
    K = malloc(ncomp * sizeof(*K));
    
    Ps = 0.0;
    Pu = 0.0;
    
    while(1) {
        status = flash_calculation_stability_analysis_QNSS(phase, K, 1e-10);
        
		/* unstable */
        if (status == 0) {
            Pu = eos->pres;
            
            if (search == 0 && !Pu_found) {
                search = 1;
			}
			else if (search == 1 && !Pu_found) {
                search = 0;
			}
                
            Pu_found = 1;
		}
        else {
            Ps = eos->pres;
            Ps_found = 1;
		}
        
        if (Pu_found && Ps_found) 
            break;
        
        if (search == 0)
            eos->pres -= search_step;
        else if (search == 1)
            eos->pres += search_step;
        
        if (eos->pres < 1.0) {
            Pu = 1.0;
            Ps = 1.0;
            break;
		}
        else if (eos->pres > P_max) {
            if (!status) {
                Pu = P_max;
                Ps = P_max;
            }
            else {
                Pu = 1.0;
                Ps = 1.0;
            }
            break;
        }
	}

    Pu_Ps[0] = Pu;
	Pu_Ps[1] = Ps;

	free(K);
    flash_calculation_phase_free(&phase);
}

/* ### Search temperature at a given pressure
# At a given pressure, search a temperature $T_s$ at which a single phase exists and a temperature $T_u$ at which two phases co-exist. The stability test is used in this function.
*/
void flash_calculation_saturation_calculation_search_Tu_Ts(EOS *eos, double *z, double P, 
        double T_est, int search, double search_step, int Ts_found, int Tu_found, 
        double *Tu_Ts)
{
	/* default values:  T_est: 100.0, search: 0 upper, 1 down, search_step: 1.0, Ts_found: 0, Tu_found: 0. */
    int ncomp = eos->ncomp;
	double Ts, Tu, *K;
	PHASE *phase;
	int status;

	eos->temp = T_est;
    eos->pres = P;
    phase = flash_calculation_phase_new(eos, z);

    K = malloc(ncomp * sizeof(*K));
    
    Ts = 0.0;
    Tu = 0.0;
    
    while(1) {
        status = flash_calculation_stability_analysis_QNSS(phase, K, 1e-10);
        
        if (status == 0) {
            Tu = eos->temp;
            
            if (search == 0 && !Tu_found) 
                search = 1;
            else if (search == 1 && !Tu_found)
                search = 0;
                
            Tu_found = 1;
		}
        else {
            Ts = eos->temp;
            Ts_found = 1;
		}
        
        if (Tu_found && Ts_found)
            break;
        
        if (search == 0) 
            eos->temp -= search_step;
        else if (search == 1)
            eos->temp += search_step;
        
		if (eos->temp < 1.0 || eos->temp > 1000.0) {
            Tu = 1.0;
            Ts = 1.0;
            break;
		}
	}

	Tu_Ts[0] = Tu;
	Tu_Ts[1] = Ts;

    free(K);
    flash_calculation_phase_free(&phase);
}

/* ### Calculate the pressure at the single-phase and two-phase boundary
# With the $P_s$ and $P_u$, the bisection method is used to search the pressure at the single-phase and two-phase boundary. The method can also be used to generate phase diagram.
*/
double flash_calculation_saturation_calculation_search_boundary_pressure(EOS *eos, double *z, double T, 
	double Pu, double Ps, int search, double tol, int search_status, double *K)
{
	/* search: 0 upper, 1 down */
	/* search_status: -1 anyone, 0 unstable, 1 stable */
	int ncomp = eos->ncomp, status;
	double P_max, P_min, P_mid;
	PHASE *phase;

    if (search == 0) {
        P_max = Ps;
        P_min = Pu;
	}
    else if (search == 1) {
        P_max = Pu;
        P_min = Ps;
	}
    else {
        printf("Saturation search boundary pressure: wrong search!");
	}

    phase = flash_calculation_phase_new(eos, z);
    
    while(1) {
        P_mid = (P_max + P_min) * 0.5;
        
        eos->pres = P_mid;
        
        status = flash_calculation_stability_analysis_QNSS(phase, K, 1e-10);
        
        if (status == 0) {
            if (search == 0) {
                P_min = P_mid;
			}
            else if (search == 1) {
                P_max = P_mid;
			}
		}
        else if (status == 1) {
            if (search == 0)
                P_max = P_mid;
            else if (search == 1)
                P_min = P_mid;
		}
        else {
            printf("Saturation search boundary pressure: wrong stability status!");
		}

        if (fabs(P_max - P_min) < tol) {
            if (search_status == -1 || search_status == status)
                break;
		}
	}

    flash_calculation_phase_free(&phase);

    return eos->pres;
}

/*  ### Calculate Y
# $$
# Y_i = z_i K_i
# $$
*/

void flash_calculation_saturation_calculation_calculate_X(double *K, double *z, double *X, int ncomp)
{
	int i;

    for (i = 0; i < ncomp; i++) {
        X[i] = z[i] * K[i];
	}
}

/* ### Calculate derivative of Y with respect to K
# $$
# \frac{\partial Y_i}{\partial K_i} = z_i
# $$
*/

void flash_calculation_aturation_calculation_calculate_X_derivative(double *z, double *dX_dK, int ncomp)
{
	int i;

    for (i = 0; i < ncomp; i++) {
        dX_dK[i] = z[i];
	}
}

/* ### Calculate fugacity-ratio corrections
# $$
# R_i = \frac{f_{z,i}}{f_{y,i}} (\sum_{j = 1}^{N_c}{Y_i})^{-1}
# $$
*/

void flash_calculation_saturation_calculation_calculate_fugacity_ratio(PHASE *phase, PHASE *phase_t, 
	double *Y, double *Ri)
{
	int i, ncomp = phase->ncomp;
	double sum_Y = 0.0;
    
    for (i = 0; i < ncomp; i++) {
        sum_Y += Y[i];
	}
    
    for (i = 0; i < ncomp; i++) {
        Ri[i] = phase->fug[i] / phase_t->fug[i] / sum_Y;
	}
}

/* ### Update incipient-phase mole numbers with the fungacity-ratio corrections
# $$
# Y_i^{n+1} = Y_i^n R_i
# $$
*/
void flash_calculation_saturation_calculation_update_X_with_fugacity_ratio(double *Y, double *Ri, int ncomp)
{
	int i;
    
    for (i = 0; i < ncomp; i++) {
        Y[i] = Y[i] * Ri[i];
	}
}

/* ### Calculate the Q value
# The recommended approach for determining saturation pressure is based on an approach proposed by Michelsen; he uses the condition:
# $$
# Q(P_{\text{sat}}, y) = 1 - \sum_{i=1}^{N_c}{z_i \frac{\phi_i(z)}{\phi_i(y)}} = 0
# $$
# Then we can have
# $$
# Q(P_{\text{sat}}, y) = 1 - \sum_{i=1}^{N_c}{y_i \frac{f_{z,i}}{f_{y,i}}} = 1 - \sum_{i=1}^{N_c}{Y_i}
# $$
*/
double flash_calculation_saturation_calculation_calculate_Q_value(double *Y, int ncomp)
{
	int i;
    double sum_Y = 0.0;
    
    for (i = 0; i < ncomp; i++) {
        sum_Y += Y[i];
	}
    
    sum_Y = 1.0 - sum_Y;
    
    return sum_Y;
}

/* ### Calculate the derivative of Q with respect to pressure
# $$
# \frac{\partial Q}{\partial P} = \sum_{i=1}^{N_c}{Y_i R_i 
#     (\frac{\partial f_{y,i}}{\partial P} \frac{1}{f_{y,i}}
#     - \frac{\partial f_{z,i}}{\partial P} \frac{1}{f_{z,i}})}
# $$
*/
double flash_calculation_saturation_calculation_calculate_Q_derivative(PHASE *phase, PHASE *phase_t, 
	double *Y, double *Ri)
{
	int i, ncomp = phase->ncomp;
    double dQ = 0.0;
    
    for (i = 0; i < ncomp; i++) {
        dQ += Y[i] * Ri[i] * (phase_t->dfug[i] / phase_t->fug[i] - phase->dfug[i] / phase->fug[i]);
	}

    return dQ;
}

/* ### Saturation Pressure Calculation
# The algorithm is from the book "Phase Behavior" by Curtis H. Whitson and Michael R. Brule.
# 1. Guess a saturation type: bubble- or dewpoint. An incorrect guess will not affect convergence, but the final K values may be "upside down".
# 2. Guess a pressure p*.
# 3. Perform Michelsen's stability test at p*.
# 4. 
#     (a) If the mixture is stable for the current value of p*,this pressure represents p* the upper bound of the search for a saturation 
		pressure on the upper curve of the phase envelope. Return to Step 1 and try a lower pressure to look for an unstable condition.
#     (b) With an unstable condition at p*, this pressure represents the lower bound in the search for a saturation pressure on the upper 
			curve of the phase envelope.
# 5. Having found an unstable solution, use the K values from the stability test to calculate incipient-phase mole numbers at bubble and 
			dewpoint with 
# $$
#                 Y_i = z_i K_i
# $$
# and 
# $$
# Y_i = z_i / K_i
# $$
# 6. Calculate the normalized incipient-phase compositions
# $$
#             y_i = Y_i / \sum{Y_i}
# $$
# 7. Calculate phase Z factors, Z_z and Z_y, and component fugacities, f_z and f_y, from the EOS at the pressure saturation-pressure estimate.
			When multiple Z-factor roots are found for a given phase, the root giving the lowest Gibbs energy should be chosen.
# 8. Calculate fugacity-ratio corrections:
# $$
#               R_i = \frac{f_{z,i}}{f_{y,i}} \sum(Y_i)^(-1)
# $$
# 9. Update incipient-phase mole numbers with the fugacity-ratio corrections:
# $$                    
#                     Y_i = Y_i * R_i^\lambda
# $$
# where four iterations use successive substitution ($\lambda = 1$) followed by a GDEM promotion with lambda given by
# $$
#                     \lambda = \| \frac{b_{11}}{b_{11} - b_{01}}\|
# $$
# where $b_{01} = \sum{\ln{R_i^n} \ln{R_i^{n-1}}}$ and $b_{11} = \sum{ln{R_i^{n-1}} \ln{R_i^{n-1}}}$.
# 10. Calculate a new estimate of saturation pressure using a Newton update:
# $$
#                     P_{\text{sat}}^{n+1} = P_{\text{sat}}^n - \frac{Q^n}{dQ/dp^n}
# $$
# If searching for an upper saturation pressure, the new pressure estimate must be higher than p*. If the new estimate is lower than p*, 
			go to Step 1 and use a new initial-pressure estimate higher than the pressure p* value.
# 11. Check for convergence. Zick suggests the following two criteria (a) 
# $$
# | 1 - \sum{Y_i} | < 10^{-13}
# $$ 
# (b) 
# $$(\sum{\ln{R_i} / \ln{Y_i/z_i}})^2 < 10^{-8}$$
# In addition, check for a trivial solution using the criterion
# $$
#                     \sum{\ln{Y_i/z_i}^2} < 10^{-4}
# $$
# 12. (a) If convergence is not achieved, return to Step 6.
#     (b) If convergence is achieved, determine the saturation type by comparing the mole fraction of the heaviest component in the mixture 
			with that in the incipient phase, where $y_N < z_N$ indicates a bubble point with $K_i = y_i / z_i$ and $y_N > z_N$ indicates a 
			dew point with $K_i = z_i / y_i$, or by comparing the density of the incipient phase with that of the feed.
*/
double flash_calculation_saturation_calculation(EOS *eos, double *z, double T, double P_est, 
        int search, double search_step, double P_max)
{
	/* search: 0 upper, 1 down */
	int i, itr, ncomp = eos->ncomp;
	double Pu_Ps[2], Pu, Ps, P0;
	double *X, *x, *Ri, *K;
    PHASE *phase_x, *phase_z;

    flash_calculation_saturation_calculation_search_Pu_Ps(eos, z, T, P_est, 
            search, search_step, 0, 0, Pu_Ps, P_max);
    Pu = Pu_Ps[0];
	Ps = Pu_Ps[1];

    if (Pu <= 1.0 && Ps <= 1.0) {
        return 1.0;
	}

    if (Pu >= P_max && Ps >= P_max) {
        return P_max;
    }
    
	/* Search unstable pressure */
    K = malloc(ncomp * sizeof(*K));
    P0 = flash_calculation_saturation_calculation_search_boundary_pressure(eos, z, T, 
            Pu, Ps, search, 0.05, 0, K);
    
    if (P0 <= 1.0 && search == 1) {
        return P0;
	}
    
    /* Initial composition */
    X = malloc(ncomp * sizeof(*X));
    x = malloc(ncomp * sizeof(*x));
    Ri = malloc(ncomp * sizeof(*Ri));

    flash_calculation_saturation_calculation_calculate_X(K, z, X, ncomp);
    
    /* Initial phase */
    phase_x = flash_calculation_phase_new(eos, x);
    phase_z = flash_calculation_phase_new(eos, z);
    eos->pres = P0;
    
    itr = 0;
    while(1) { 
		double Q, dQ, dp;

        flash_calculation_calculate_trial_phase_composition(X, x, ncomp);
    
        /* Calculate phase z fugacity */
        flash_calculation_compute_phase_parameter(phase_z);
        flash_calculation_calculate_compressibility_factor(phase_z);
        flash_calculation_calculate_fugacity(phase_z);
        
        /* Calculate phase x fugacity */
        flash_calculation_compute_phase_parameter(phase_x);
        flash_calculation_calculate_compressibility_factor(phase_x);
        flash_calculation_calculate_fugacity(phase_x);
        
        /* Calcualte fugacity-ratio corrections Ri */
        flash_calculation_saturation_calculation_calculate_fugacity_ratio(phase_z, phase_x, 
                X, Ri);
        
		/* # Update incipient-phase mole numbers with
           # fugacity-ratio corrections */
        flash_calculation_saturation_calculation_update_X_with_fugacity_ratio(X, Ri, ncomp);
        
        /* # Calculate a new estimate of saturation pressure
           # using a Newton-Raphson update */
        Q = flash_calculation_saturation_calculation_calculate_Q_value(X, ncomp);
        dQ = flash_calculation_saturation_calculation_calculate_Q_derivative(phase_z, phase_x, X, Ri);
        dp = - Q / dQ;
        eos->pres += dp;

        if (eos->pres < 1.0) {
            eos->pres -= dp;

            eos->pres = (eos->pres + 1.0) / 2.0;
        }
    
        /* Check convergence */
        if (fabs(Q) < 1e-10) 
            break;
            
        itr += 1;
        if (itr > 1000) {
            //printf("##### WARNING: Saturation calculation reach maximum iterations!\n");
            break;
		}
	}
      
    free(K);
    free(X);
    free(x);
    free(Ri);

    flash_calculation_phase_free(&phase_x);
    flash_calculation_phase_free(&phase_z);

    return eos->pres;
}

