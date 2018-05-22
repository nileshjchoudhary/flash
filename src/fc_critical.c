#include "fc.h"
/* # ### Calculate the derivative of Wilson's equation with respect to pressure */
/*
   K[i] = EOS.comp[i].PC / P \
 * np.exp(5.37 * (1.0 + EOS.comp[i].AC) \
 * (1.0 - EOS.comp[i].TC / T))
 dK[i]/dp = - Pc / (P * P) * Exp()
 */      
double flash_calculation_calculate_initial_K_derivative_with_pressure(EOS *eos)
{
    int i, ncomp = eos->ncomp;
    double P, T, dK;

    P = eos->pres;
    T = eos->temp;

    dK = 0.0;
    for (i = 0; i < ncomp; i++) {
        double temp;

        temp = exp(5.37 * (1.0 + eos->comp_list->comp[i].AC) * (1.0 - eos->comp_list->comp[i].TC / T));
        dK += - eos->comp_list->comp[i].PC / (P * P) * temp;
    }

    return dK;
}

/* ### Calculate the derivative of Wilson's equation with respect to temperature */
/*
   dK[i]/dT = Pc / P * Exp() * (5.37 * (1.0 + Ac)) * Tc / (T * T)
   */
double flash_calculation_calculate_initial_K_derivative_with_temperature(EOS *eos)
{
    int i, ncomp = eos->ncomp;
    double P, T, dK;
    COMP *comp = eos->comp_list->comp;

    P = eos->pres;
    T = eos->temp;

    dK = 0.0;
    for (i = 0; i < ncomp; i++) {
        double temp;

        temp = exp(5.37 * (1.0 + comp[i].AC) * (1.0 - comp[i].TC / T));
        dK += comp[i].PC / P * temp * 5.37 * (1.0 + comp[i].AC)
            * comp[i].TC / (T * T);
    }

    return dK;
}

/* ### Calculate the pressure at given compositions, temperature and $F_v$
# The equilibirum equation and Richford-Rice equation are used in this calculation.
*/
/*
   Equilibrium equation: 
   G[i] = ln(K[i]) + ln(phi_v[i]) - ln(phi_l[i]) = 0, 
   for i = 0, ..., Nc
   Material balance:
   G[Nc+1] = Sum_i((K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0)))
   where Ki = x_v[i] / x_l[i]
   with primary variables P and K[i].
   */
double * flash_calculation_calculate_status_at_fixed_temperature_and_Fv(EOS *eos, double *z, 
        double T, double Fv, double P_max, double P_min, double P_est, double *K)
{
    /* default values:   P_max: 1000.0, P_min: 1.0, P_est: -1.0, K: NULL */
    int itr, ncomp = eos->ncomp, i;
    double F, *x_l, *x_v, *G, *dG, error, *dv, dK;
    PHASE *phase_L, *phase_V;
    int break_flag;

    /* Initial guess */
    if (K == NULL && P_est < 0.0) {
        double Pu_Ps[2];

        K = malloc(ncomp * sizeof(*K));

        flash_calculation_saturation_calculation_search_Pu_Ps(eos, z, T, 50.0, 0, 1.0, 1, 0, Pu_Ps);
        eos->pres = Pu_Ps[0];
        eos->temp = T;
        flash_calculation_estimate_K(eos, K);
        itr = 0;
        while(1) {
            F = flash_calculation_calculate_RachfordRice_equation_value(K, z, Fv, ncomp);

            if (fabs(F) < 1e-3)
                break;

            dK = flash_calculation_calculate_initial_K_derivative_with_pressure(eos);
            eos->pres += - F / dK;

            flash_calculation_estimate_K(eos, K);

            itr += 1;
            if (itr > 30) {
                break;
            }
        }

        if (eos->pres < 0.0) {
            eos->pres = Pu_Ps[0];
            flash_calculation_estimate_K(eos, K);
        }
    }
    else if (K == NULL && P_est > 0.0) {
        K = malloc(ncomp * sizeof(*K));

        eos->pres = P_est;
        eos->temp = T;
        flash_calculation_estimate_K(eos, K);
    }
    else {
        eos->pres = P_est;
        eos->temp = T;
    }

    /* Initial compositions x_v and x_l */
    x_l = malloc(ncomp * sizeof(*x_l));
    x_v = malloc(ncomp * sizeof(*x_v));

    /* Initial liquid and vapour phase */
    phase_L = flash_calculation_phase_new(eos, x_l);
    phase_V = flash_calculation_phase_new(eos, x_v);

    /* Initial G */
    G = malloc((ncomp + 1) * sizeof(*G));
    dG = malloc((ncomp + 1) * (ncomp + 1) * sizeof(*dG));

    /* Calculate liquid and vapour composition */
    flash_calculation_calculate_composition(K, z, Fv, x_l, x_v, ncomp);

    itr = 0;
    break_flag = 0;
    while(1) {
        /* Calculate phase vapour */
        flash_calculation_compute_phase_parameter(phase_V);
        flash_calculation_calculate_compressibility_factor(phase_V);
        flash_calculation_calculate_fugacity(phase_V);

        /* Calculate phase liquid */
        flash_calculation_compute_phase_parameter(phase_L);
        flash_calculation_calculate_compressibility_factor(phase_L);
        flash_calculation_calculate_fugacity(phase_L);

        /* Calculate the equations */
        error = flash_calculation_phase_diagram_construction_calculate_equations(phase_L,
                phase_V, K, z, Fv, G);

        /* Check if convergence */
        if (error < 1e-10) 
            break;

        /* Calculate the derivative of the equations */
        flash_calculation_phase_diagram_construction_calculate_equations_derivative(phase_L,
                phase_V, K, z, Fv, dG);

        /* Update variables */
        dv = flash_calculation_phase_diagram_construction_update_variables(dG, G, ncomp + 1);
        for (i = 0; i < ncomp; i++) {
            K[i] += dv[i];
        }
        eos->pres += dv[ncomp];

        flash_calculation_calculate_composition(K, z, Fv, x_l, x_v, ncomp);

        if (eos->pres > P_max) {
            if (break_flag) {
                eos->pres = -1.0;
                free(dv);
                break;
			}
            else {
                break_flag = 1;
                eos->pres -= dv[ncomp];
                eos->pres = (P_max + eos->pres) * 0.5;
			}
		}

        if (eos->pres < P_min) {
            if (break_flag) {
                eos->pres = -1.0;
                free(dv);
                break;
			}
            else {
                break_flag = 1;
                eos->pres -= dv[ncomp];
                eos->pres = (P_min + eos->pres) * 0.5;
			}
		}
                
		free(dv);

        itr += 1;
        if (itr > 200) {
            //printf("WARNING: Calculate_status_at_fixed_temperature_and_Fv reach maximum iterations!\n");
            if (error > 1e-4) {
                eos->pres = -1.0;
			}
			break;
		}
	}

    free(x_l);
    free(x_v);
    free(G);
    free(dG);
    flash_calculation_phase_free(&phase_L);
    flash_calculation_phase_free(&phase_V);
    
    return K;
}


/* ### Calculate the temperature at given compositions, pressure and $F_v$
# The equilibirum equation and Richford-Rice equation are used in this calculation.
*/
/*
        Equilibrium equation: 
            G[i] = ln(K[i]) + ln(phi_v[i]) - ln(phi_l[i]) = 0, 
                for i = 0, ..., Nc
        Material balance:
            G[Nc+1] = Sum_i((K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0)))
        where Ki = x_v[i] / x_l[i]
        with primary variables T and K[i].
*/

double * flash_calculation_calculate_status_at_fixed_pressure_and_Fv(EOS *eos, double *z, 
        double P, double Fv, double T_max, double T_min, double T_est, double *K) 
{
    /* default values:   T_max: 1000.0, T_min: 1.0, T_est: -1.0, K: NULL */
    int ncomp = eos->ncomp, itr;
	double Tu_Ts[2], F, dK, *x_l, *x_v, *G, *dG, error, *dv;
	PHASE *phase_L, *phase_V;
	int break_flag, i;
    
    /* Initial guess */
    if (K == NULL && T_est < 0.0) {
		K = malloc(ncomp * sizeof(*K));

        flash_calculation_saturation_calculation_search_Tu_Ts(eos, z, P,  50.0, 0, 1.0, 1, 0, Tu_Ts);
        eos->pres = P;
        eos->temp = Tu_Ts[0];
        flash_calculation_estimate_K(eos, K);
        itr = 0;

        while(1) {
            F = flash_calculation_calculate_RachfordRice_equation_value(K, z, Fv, ncomp);
            
            if (fabs(F) < 1e-3)
                break;
                
            dK = flash_calculation_calculate_initial_K_derivative_with_temperature(eos);
            eos->temp += - F / dK;
            
            flash_calculation_estimate_K(eos, K);
            
            itr += 1;
            if (itr > 30)
                break;
		}

        if (eos->temp < 0.0) {
            eos->temp = Tu_Ts[0];
            flash_calculation_estimate_K(eos, K);
		}
	}
    else if (K == NULL && T_est > 0.0) {
		K = malloc(ncomp * sizeof(*K));

        eos->pres = P;
        eos->temp = T_est;
        flash_calculation_estimate_K(eos, K);
	}
    else {
        eos->pres = P;
        eos->temp = T_est;
	}

    /* Initial compositions x_v and x_l */
    x_l = malloc(ncomp * sizeof(*x_l));
    x_v = malloc(ncomp * sizeof(*x_v));
    
    /* Initial liquid and vapour phase */
    phase_L = flash_calculation_phase_new(eos, x_l);
    phase_V = flash_calculation_phase_new(eos, x_v);
    
    /* Initial G */
    G = malloc((ncomp + 1) * sizeof(*G));
    dG = malloc((ncomp + 1) * (ncomp + 1) * sizeof(*dG));
    
    /* Calculate liquid and vapour composition */
    flash_calculation_calculate_composition(K, z, Fv, x_l, x_v, ncomp);
    
    itr = 0;
    break_flag = 0;
    while(1) {
        /* Calculate phase vapour */
        flash_calculation_compute_phase_parameter(phase_V);
        flash_calculation_calculate_compressibility_factor(phase_V);
        flash_calculation_calculate_fugacity(phase_V);
        
        /* Calculate phase liquid */
        flash_calculation_compute_phase_parameter(phase_L);
        flash_calculation_calculate_compressibility_factor(phase_L);
        flash_calculation_calculate_fugacity(phase_L);
        
        /* Calculate the equations */ 
        error = flash_calculation_phase_diagram_construction_calculate_equations(phase_L,
			phase_V, K, z, Fv, G);
            
        /* Check if convergence */
        if (error < 1e-10)
            break;

        /* Calculate the derivative of the equations */
        flash_calculation_phase_diagram_construction_calculate_equations_derivative_with_K_T(phase_L, 
			phase_V, K, z, Fv, dG);
            
        /* Update variables */
        dv = flash_calculation_phase_diagram_construction_update_variables(dG, G, ncomp + 1);
        for (i = 0; i < ncomp; i++) 
            K[i] += dv[i];
        eos->temp += dv[ncomp];
        
        flash_calculation_calculate_composition(K, z, Fv, x_l, x_v, ncomp);
        
        if (eos->temp > T_max) {
            if (break_flag) {
                eos->temp = -1.0;
                free(dv);
                break;
			}
            else {
                break_flag = 1;
                eos->temp -= dv[ncomp];
                eos->temp = (T_max + eos->temp) * 0.5;
			}
		}
        
        if (eos->temp < T_min) {
            if (break_flag) {
                eos->temp = -1.0;
                free(dv);
                break;
			}
            else {
                break_flag = 1;
                eos->temp -= dv[ncomp];
                eos->temp = (T_min + eos->temp) * 0.5;
			}
		}

        free(dv);
        itr += 1;
        if (itr > 100)
            break;
	}

    free(x_l);
    free(x_v);
    free(G);
    free(dG);

    flash_calculation_phase_free(&phase_L);
    flash_calculation_phase_free(&phase_V);
    
    return K;	
}

/* ### Initial temperature guess for critical point calcuation
# $$
# T_{\text{init}} = \sum{z_i T_{c,i}}
# $$
*/
 /*  """
        The initial guess Ti = 1.5 * sum(x_i * Tc_i)
    """
 */

double flash_calculation_critial_point_calculate_init_T_guess(EOS *eos, double *z)
{
    int i, ncomp = eos->ncomp;
	double Ti;
	COMP *comp = eos->comp_list->comp;
    
    Ti = 0.0;
    for (i = 0; i < ncomp; i++) {
        Ti += z[i] * comp[i].TC;
	}
    
    Ti *= 1.3;
    
    return Ti;
}

/* ### Calculate the volume functions $F_1(\kappa)$ -- $F_8(\kappa)$, $\kappa = \frac{v}{b}$
# $$
# F_1 = \frac{1}{\kappa - 1} \\
# F_2 = 2 \frac{\frac{\sigma_1}{\kappa + \sigma_1} - \frac{\sigma_2}{\kappa + \sigma_2}}{\sigma_1 - \sigma_2}\\
# F_3 = \frac{(\frac{\sigma_1}{\kappa + \sigma_1})^2 - (\frac{\sigma_2}{\kappa + \sigma_2})^2}{\sigma_1 - \sigma_2} \\
# F_4 = \frac{(\frac{\sigma_1}{\kappa + \sigma_1})^3 - (\frac{\sigma_2}{\kappa + \sigma_2})^3}{\sigma_1 - \sigma_2} \\
# F_5 = 2 \frac{\ln{\frac{\kappa + \sigma_1}{\kappa + \sigma_2}}}{\sigma_1 - \sigma_2}\\
# F_6 = F_2 - F_5 \\
# F_7 = - \frac{F_2}{1 + F_1} \\
# F_8 = \frac{F_3}{1 + F_1} 
# $$
*/
void flash_calculation_critical_point_volume_functions(double kappa, double sigma1, double sigma2, 
	double *F)
{
	double temp1, temp2, down;

    F[0] = 1.0 / (kappa - 1.0);
    
    temp1 = sigma1 / (kappa + sigma1);
    temp2 = sigma2 / (kappa + sigma2);
    down = sigma1 - sigma2;
    
    F[1] = 2.0 * (temp1 - temp2) / down;
    F[2] = (temp1 * temp1 - temp2 * temp2) / down;
    F[3] = (temp1 * temp1 * temp1 - temp2 * temp2 * temp2) / down;
    F[4] = 2.0 * log((kappa + sigma1) / (kappa + sigma2)) / down;
    F[5] = F[1] - F[4];
    F[6] = - F[1] / (1.0 + F[0]);
    F[7] = F[2] / (1.0 + F[0]);
}

/* ### Calculate $\alpha$ and $\beta$
# $$
# \beta_i = \frac{b_i}{b}\\
# \alpha_i = \frac{\sum_j{y_j a_{ij}}}{a}
# $$
*/

void flash_calculation_critical_point_calculate_alpha_beta(PHASE *phase, double *alpha, double *beta)
{
    int i, j, ncomp = phase->ncomp;
	COMP *comp = phase->eos->comp_list->comp;
    
    for (i = 0; i < ncomp; i++) {
		double aij;

		aij = 0.0;
        alpha[i] = 0.0;
        
        for (j = 0; j < ncomp; j++) {
            aij = sqrt(phase->ai[i] * phase->ai[j]) * (1.0 - comp[i].binary[j]);
            alpha[i] += phase->mf[j] * aij;
		}
        alpha[i] /= phase->a;
        
        beta[i] = phase->bi[i] / phase->b;
	}
}


/* ### Calculate the Q matrix
# $$
# Q_{ij} = R T (\frac{\partial \ln{f_i}}{\partial N_j})
# $$
# and 
# $$
# RT\frac{\partial \ln{f_i}}{\partial N_j} = RT[\frac{\delta_{ij}}{y_i} + (\beta_i + \beta_j) F_1 + \beta_i \beta_j F_1^2]
#     + \frac{a}{b}[\beta_i \beta_j F_3 - a_{ij}/a F_5 + (\beta_i \beta_j - \alpha_i \beta_j - \alpha_j \beta_i) F_6]
# $$
*/

void flash_calculation_critical_point_calculate_Q_matrix(double *F, PHASE *phase, double *alpha, double *beta, double *Q)
{
    int i, j, ncomp = phase->ncomp;
	double sigmaij, Qij, aij;
	EOS *eos = phase->eos;
	COMP *comp = phase->eos->comp_list->comp;
    
    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            sigmaij = 0;
            
            if (i == j)
                sigmaij = 1.0;
            else
                sigmaij = 0.0;
                
            aij = sqrt(phase->ai[i] * phase->ai[j])
				* (1.0 - comp[i].binary[j]);
                
            Qij = phase->R * eos->temp
				* (sigmaij / phase->mf[i] + (beta[i] + beta[j]) * F[0]
					+ beta[i] * beta[j] * F[0] * F[0]);
                
            Qij += phase->a / phase->b * (beta[i] * beta[j] * F[2]
				- aij / phase->a * F[4] 
				+ (beta[i] * beta[j] - alpha[i] * beta[j] 
				- alpha[j] * beta[i]) * F[5]);
            
			//Q[i * ncomp + j] = Qij;
			Q[i * ncomp + j] = Qij / phase->R / eos->temp;
		}
	}
}

double flash_calculation_calculate_matrix_determinant(double *Q, int dim)
{
    int i, j, k, index;
    double ai1, signal = 1.0, minus = -1.0, *Qi;
    double det, deti;

    if (dim == 1) {
        return *Q;
    }

    Qi = malloc((dim - 1) * (dim - 1) * sizeof(*Qi));

    det = 0.0;
    for (i = 0; i < dim; i++) {
        ai1 = Q[i * dim + 0];

        index = 0;
        for (j = 0; j < dim; j++) {
            if (i == j) 
                continue;

            for (k = 0; k < dim - 1; k++) {
                Qi[index * (dim - 1) + k] = Q[j * dim + k + 1];
            }
            index += 1;
        }

        deti = flash_calculation_calculate_matrix_determinant(Qi, dim - 1);
        det += ai1 * signal * deti;

        signal *= minus;
    }

    free(Qi);

    return det;
}

/* ### Calculate the determinant of matrix Q */
/*
    """
        Calculate the matrix Q and the determinant
    """
*/
double flash_calculation_critical_point_calculate_Q_det(double *F, PHASE *phase, 
        double Tc, double *Qij)
{   
	int i, j, ncomp = phase->ncomp;
	double *alpha, *beta, det_Q;
    
    alpha = malloc(ncomp * sizeof(*alpha));
    beta = malloc(ncomp * sizeof(*beta));

    phase->eos->temp = Tc;
    
    /* Calculate ai, bi, a, b */
    flash_calculation_compute_phase_parameter(phase);
    
    /* Calculate alpha and beta */
    flash_calculation_critical_point_calculate_alpha_beta(phase, alpha, beta);
    
    /* Calculate Q matrix */
    flash_calculation_critical_point_calculate_Q_matrix(F, phase, alpha, beta, Qij);
    
	det_Q = flash_calculation_calculate_matrix_determinant(Qij, ncomp);

    free(alpha);
    free(beta);
    
    return det_Q;
}

/* ### Calculate the derivative of Q determinant with respect to T
# $$
# \frac{\partial |Q|}{\partial T} = \frac{|Q(T + \delta T)| - |Q(T - \delta T)|}{2 \delta T}
# $$
*/
/*
    """
        Numerical method is used to calculate the direvative:
            dQ/dT = (Q_(i+1) - Q_(i-1)) / (2.0 * dT),
        where dT = Tc * 1e-7, Q_(i+1) = Q(Tc + dT), and 
        Q_(i-1) = Q(Tc - dT).
    """
*/
double flash_calculation_critical_point_calculate_Q_det_derivative(double *F, 
        PHASE *phase, double Tc)
{   
    int ncomp = phase->ncomp;
	double dT = Tc * 1e-9, det_Q_1, det_Q_2, d_det_Q;
    double *Q1, *Q2;
    
    Q1 = malloc(ncomp * ncomp * sizeof(*Q1));
    Q2 = malloc(ncomp * ncomp * sizeof(*Q2));

    det_Q_1 = flash_calculation_critical_point_calculate_Q_det(F, phase, Tc - dT, Q1);
    det_Q_2 = flash_calculation_critical_point_calculate_Q_det(F, phase, Tc + dT, Q2);
    
    d_det_Q = (det_Q_2 - det_Q_1) / (2.0 * dT);

    free(Q1);
    free(Q2);
    
    return d_det_Q;
}

/* ### Calculate the cubic form
# $$
# C = RT[- \sum_i{\frac{\Delta N_i^3}{y_i^2}} - 3\bar{N}(\bar{\beta} F_1)^2 + 2 (\bar{\beta}F_1)^3]
#      + \frac{a}{b} [3 \bar{\beta}^2 (2 \bar{\alpha} - \bar{\beta}) (F_3 + F_6) 
#                      - 2 \bar{\beta}^3 F_4 - 3 \bar{\beta} \bar{a} F_6]
# $$
# where
# $$
# \bar{N} = \sum_i{\Delta N_i} \\
# \bar{\beta} = \sum_i{\Delta N_i \beta_i} \\
# \bar{\alpha} = \sum_i{\Delta N_i \alpha_i} \\
# \bar{a} = \frac{1}{a} \sum_i{\sum_j{\Delta N_i \Delta N_j a_{ij}}}
# $$
*/
double flash_calculation_critical_point_calculate_cubic_form(PHASE *phase, 
        double *N, double *alpha, double *beta, double kappa)
{
    int i, j, ncomp = phase->ncomp;
    double *F, N_b, beta_b, alpha_b, a_b, aij, C;
    EOS *eos = phase->eos;
    COMP *comp = eos->comp_list->comp;
    
    F = malloc(8 * sizeof(*F));
    flash_calculation_critical_point_volume_functions(kappa, eos->para_sigma1, 
            eos->para_sigma2, F);
    
    /* Calculate the parameters for the cubic form calculation */
    N_b = 0.0;
    beta_b = 0.0;
    alpha_b = 0.0;
    a_b = 0.0;

    for (i = 0; i < ncomp; i++) {
        N_b += N[i];
        beta_b += N[i] * beta[i];
        alpha_b += N[i] * alpha[i];
        
        for (j = 0; j < ncomp; j++) {
            aij = sqrt(phase->ai[i] * phase->ai[j]) 
                * (1.0 - comp[i].binary[j]);
            a_b += aij * N[i] * N[j];
        }
    }
    a_b = a_b / phase->a;
    
    C = 0.0;
    for (i = 0; i < ncomp; i++) {
        C += - pow(N[i], 3.0) / pow(phase->mf[i], 2.0);
    }
        
    C += 3.0 * N_b * pow(beta_b * F[0], 2.0)
        + 2.0 * pow(beta_b * F[0], 3.0);
    C += phase->a / (phase->b * phase->R * eos->temp)
        * (3.0 * beta_b * beta_b * (2.0 * alpha_b - beta_b) * (F[2] + F[5])
                - 2.0 * pow(beta_b, 3.0) * F[3] - 3.0 * beta_b * a_b * F[5]);

    free(F);
        
    //return C * phase->R * eos->temp;
    return C;
}


/* ### Calculate the derivative of the cubic form C with respect to $\kappa$
# $$
# \frac{\partial C}{\partial \kappa} = \frac{C(\kappa + \delta \kappa) - C(\kappa - \delta \kappa)}{2 \delta \kappa}
# $$
*/

double flash_calculation_critical_point_calculate_cubic_form_derivative(PHASE *phase, 
        double *N, double *alpha, double *beta, double kappa)
{
    double dC, C1, C2, dkappa;

    dkappa = kappa * 1e-9;
    
    C1 = flash_calculation_critical_point_calculate_cubic_form(phase, N, alpha, 
            beta, kappa - dkappa);
    C2 = flash_calculation_critical_point_calculate_cubic_form(phase, N, alpha, 
            beta, kappa + dkappa);
    
    dC = (C2 - C1) / (2.0 * dkappa);
    
    return dC;
}

/* 
        ### From the following equation
        ###     (Qr  QN)(xr)   (0)
        ###     (QN' Q0)(xN) = (0),
        ### we have Qr * xr + QN * xN = 0
        ###     and QN' * xr + Q0 * xN = 0.
        ### To solve xr, we need to solve a linear system
        ###         Qr * xr = - QN * xN
        ### Let xN = 1.0 and solve xr, then normalize the vector x = (xr, xN) by 
        ### dividing ||x||_2
*/

void flash_calculation_critical_point_calculate_normal_N(double *Q, int ncomp, 
        double *x)
{
    int i, j;
    double *Qr, *QN, *b, *xr, sum, x_N, x_norm;

    Qr = malloc((ncomp - 1) * (ncomp - 1) * sizeof(*Qr));
    QN = malloc((ncomp - 1) * sizeof(*QN));
    b = malloc((ncomp - 1) * sizeof(*b));
    xr = malloc((ncomp - 1) * sizeof(*xr));

    for (i = 0; i < ncomp - 1; i++) {
        for (j = 0; j < ncomp - 1; j++) {
            Qr[i * (ncomp - 1) + j] = Q[i * ncomp + j];
        }

        QN[i] = Q[i * ncomp + ncomp - 1];
    }
   
    /* Set x_N = 1.0 */
    x_N = 1.0;

    /* Calculate RHS = - QN * xN  */
    for (i = 0; i < ncomp - 1; i++) {
        b[i] = - x_N * QN[i];
    }
    
    /* Solve the linear equation Qr * x = b  */
    flash_calculation_solve_dense_linear_system(Qr, b, xr, ncomp - 1);
    
    sum = 0.0;
    for (i = 0; i < ncomp - 1; i++) {
        x[i] = xr[i];
        sum += xr[i] * xr[i];
    }
    x[ncomp - 1] = x_N;
    sum += x_N * x_N;

    x_norm = sqrt(sum);
    for (i = 0; i < ncomp; i++) {
        x[i] = x[i] / x_norm;
    }

    free(Qr);
    free(QN);
    free(b);
    free(xr);
}

/* ### Critical Point Calculation
# The critical point calculation algorithm used in this code is based on the 
    paper "The Calculation of Critical Points" by Robert A. Heidemann and Ahmed M. Khalil 
    and "Calculation of Critical Points from Cubic Two-Constant Equations of State" by 
    Michael L. Michelsen and Robert A. Heidemann.
# 
# A necessary condition for a point to lie on the limit of stability is that the matrix Q with elements
# $$
# Q_{ij} = (\frac{\partial^2 A}{\partial n_j \partial n_i})
# $$
# should have a zero determinant
# $$
# Q = Det(Q(Q_{ij})) = 0
# $$
# Or equivalently, there should be a vector 
# $$
#             N = (n_1, \cdots, n_{N_c})^T
# $$
# which satisfies the equations
# $$
#             Q \cdot N = 0 \\
# C = \sum_k \sum_j \sum_i (\frac{\partial^3 A}{\partial n_k \partial n_j \partial n_i}) n_i n_j n_k = 0
# $$
# 
# Procedure:
# 
# 1. Evaluating the stability limit: $|Q| = 0$. To converge to the correct temperature, it is necessary to 
    make the initial guess high enough. The guess we use is 
# $$
# T_{\text{init}} = 1.5 \sum(x_i T_{c,i})
# $$
# At each volume $\kappa$, the temperature is found by the Newton procedure, with numerical differentiation 
    to obtain $\frac{\partial Q}{\partial T}$. In the numerical differentiation we take
# $$                
#                 \frac{\delta T}{T} = 10^{-7}
# $$
# and the criterion of convergence is that between successive iterations
# $$
# \frac{|\delta T|}{T} <= 10^{-4}
# $$
# 
# 2. Evaluation of $\Delta N$. We first take $\Delta N_n = 1$, then solve the linear system to 
    find $\Delta N_1, \cdots, \Delta N_{n-1}$. It is proved to be important to scale $\Delta N$ by dividing by 
# $$
# [\sum(\Delta N_i)]^{\frac{1}{2}}
# $$
# 
# 3. Evaluation of the Cubic Form.
# The Newton-Raphson procedure converges monotonically to the critical volume from an inital guess of 
# $$
# v = 4b,
# $$
# which means $\kappa = 4.0$.
# The numerical differentiation is also used to obtain the derivative of C with respect to $\kappa$, and we set 
# $$
# \frac{\delta \kappa}{\kappa} = 10^{-7}.
# $$
# 
# 4. The step 1-3 will be repeated until T and $\kappa$ converge.
# 
*/

CRITICAL_POINT * flash_calculation_critical_point_calculation(EOS *eos, 
        double *z, double kappa_init, double Tc_init)
{
    /* default value: kappa_init 3.5, Tc_init -1 */
    int ncomp = eos->ncomp, itr, itr_Q, itr_C;
    double *alpha, *beta, *F;
    PHASE *phase;
    double Tc, Pc, kappa, Tc0, kappa0, det_Q, d_det_Q, dTc, *Q, *x;
    double C, dC, dkappa, dTc0, dkappa0, v;
    CRITICAL_POINT *cp;

    cp = malloc(sizeof(*cp));
    
    alpha = malloc(ncomp * sizeof(*alpha));
    beta = malloc(ncomp * sizeof(*beta));
    F = malloc(8 * sizeof(*F));
    Q = malloc(ncomp * ncomp * sizeof(*Q));
    x = malloc(ncomp * sizeof(*x));
    
    /* initial phase */
    phase = flash_calculation_phase_new(eos, z);
    
    /* Initial temperature guess */
    if (Tc_init < 0.0) {
        Tc_init = flash_calculation_critial_point_calculate_init_T_guess(eos, z);
    }
    
    Tc = Tc_init;
    
    /* Initial volume functions */
    kappa = kappa_init;
    
    itr = 0;
    while(1) {
        /* Calculate volume functions */
        flash_calculation_critical_point_volume_functions(kappa, eos->para_sigma1,
                eos->para_sigma2, F);
        
        itr_Q = 0;
        Tc0 = Tc;
        /* Evaluating the stability limit */
        while(1) {
            /* Calculae determinant of Q */
            det_Q = flash_calculation_critical_point_calculate_Q_det(F, phase, Tc, Q);
            
            /* Calculate derivative of determinant of Q */
            d_det_Q = flash_calculation_critical_point_calculate_Q_det_derivative(F, phase, Tc);

            /* Update Tc */ 
            dTc = - det_Q / d_det_Q;
            Tc += dTc;

            //printf("det_Q: %e, d_det_Q: %e, dTc: %e\n", det_Q, d_det_Q, dTc);
            
            if (Tc < 0.0) {
                Tc -= dTc;
                Tc *= 0.5;
            }
        
            /* Check if converged */
            if (fabs(dTc / Tc) < 1e-5) 
                break;
            
            itr_Q += 1;
            if (itr_Q > 100)
                break;
        }
        
        eos->temp = Tc;
        
        det_Q = flash_calculation_calculate_matrix_determinant(Q, ncomp);

        /* Evaluation of \Delta N */ 
        flash_calculation_critical_point_calculate_normal_N(Q, ncomp, x);
    
        /*## Evaluation of the Cubic Form
        ## Calculate ai, bi, a, b */
        flash_calculation_compute_phase_parameter(phase);
                   
        /* Calculate alpha and beta */
        flash_calculation_critical_point_calculate_alpha_beta(phase, alpha, beta);
        itr_C = 0;
        kappa0 = kappa;
        C = 0.0;

        while(1) {
            /* Calculate the cubic form */
            C = flash_calculation_critical_point_calculate_cubic_form(phase, x, alpha,
                    beta, kappa);
        
            /* Calculate the derivative of the cubic form */
            dC = flash_calculation_critical_point_calculate_cubic_form_derivative(phase, 
                    x, alpha, beta, kappa);
        
            /* Update kappa */
            if (fabs(dC) > 1e-10) {
                dkappa = - C / dC;
            }
            else {
                break;
            }

            kappa += dkappa;

            //printf("itr[%d]  C: %e, dC: %e, dkappa: %e\n", itr, C, dC, dkappa);
            
            if (kappa < 0.0) {
                kappa -= dkappa;
                kappa *= 0.5;
            }
        
            /* Check if converged */
            if (fabs(dkappa / kappa) < 1e-5) 
                break;
        
            itr_C += 1;
            if (itr_C > 100) 
                break;
        }
        
        /* print("kappa: %e" %kappa) */
        dTc0 = fabs(Tc - Tc0);
        dkappa0 = fabs(kappa - kappa0);

#if 0
        printf("dTc0: %e, dkappa0: %e\n", dTc0, dkappa0);
        printf("Tc: %e, kappa: %e\n", Tc, kappa);
        printf("%e, %e\n", dTc0 / Tc, dkappa0 / kappa);
#endif
        
        if (((dTc0 / Tc) < 1e-4) && ((dkappa0 / kappa) < 1e-4)) 
            break;
        
        itr += 1;
        if (itr > 200) {
            //printf("#### WARNING: Critical point calculation reach maximum iterations!\n");
            break;
        }
            
        if (kappa <= 1.0) {
            if (kappa_init < 1.0) {
                Pc = -1.0;
                Tc = -1.0;
                cp->Pc = Pc;
                cp->Tc = Tc;

                return cp;
            }
                
            kappa = kappa_init - 0.1;
            kappa_init -= 0.1;
            Tc = Tc_init * 0.9;
        }
    }
        
    /* kappa = v / b */
    eos->temp = Tc;
    flash_calculation_compute_phase_parameter(phase);
    v = kappa * phase->b;
    Pc = phase->R * Tc / (v - phase->b) - phase->a 
        / ((v + eos->para_sigma1 * phase->b)
                * (v + eos->para_sigma2 * phase->b));
    
    cp->Tc = Tc;
    cp->Pc = Pc;

    free(alpha);
    free(beta);
    free(F);
    free(Q);
    free(x);

    flash_calculation_phase_free(&phase);
            
    return cp;
}

void flash_calculation_critical_point_free(CRITICAL_POINT **cp)
{
    free(*cp);
}
