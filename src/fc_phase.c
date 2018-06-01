#include "fc.h"

PHASE * flash_calculation_phase_new(EOS *eos, double *mf)
{
    PHASE *phase;
    int ncomp;

    phase = malloc(sizeof*(phase));

    ncomp = eos->ncomp;
    phase->ncomp = ncomp;
    phase->eos = eos;
    phase->mf = mf;
    phase->R = 82.05736;
    phase->density = 0.0;

    phase->A = 0.0;
    phase->dAp = 0.0;
    phase->dA_dT = 0.0;
    phase->dA_dT2 = 0.0;
    phase->dAx = malloc(ncomp * sizeof(*(phase->dAx)));

    phase->B = 0.0;
    phase->dBp = 0.0;
    phase->dB_dT = 0.0;
    phase->dB_dT2 = 0.0;
    phase->dBx = malloc(ncomp * sizeof(*(phase->dAx)));

    phase->nroot = 0;
    phase->Z = -1;
    phase->dZ = 0;
    phase->dZ_dT = 0.0;
    phase->dZ_dT2 = 0.0;
    phase->dZ_dx = malloc(ncomp * sizeof(*(phase->dZ_dx)));
    phase->phase_no = -1;

    phase->fug = malloc(ncomp * sizeof(*(phase->fug)));
    phase->phi = malloc(ncomp * sizeof(*(phase->phi)));
    phase->dphi = malloc(ncomp * sizeof(*(phase->dphi)));
    phase->dphi_dT = malloc(ncomp * sizeof(*(phase->dphi_dT)));
    phase->dphi_dx = malloc(ncomp * ncomp * sizeof(*(phase->dphi_dx)));
    phase->dfug = malloc(ncomp * sizeof(*(phase->dfug)));

    phase->ai = malloc(ncomp * sizeof(*(phase->ai)));
    phase->dai_dT = malloc(ncomp * sizeof(*(phase->dai_dT)));
    phase->dai_dT2 = malloc(ncomp * sizeof(*(phase->dai_dT2)));
    phase->a = 0.0;
    phase->da_dT = 0.0;
    phase->da = malloc(ncomp * sizeof(*(phase->da)));
    phase->dda_dT = malloc(ncomp * sizeof(*(phase->dda_dT)));
    phase->dda_dT2 = malloc(ncomp * sizeof(*(phase->dda_dT2)));

    phase->bi = malloc(ncomp * sizeof(*(phase->dZ_dx)));
    phase->b = 0.0;
    phase->db = malloc(ncomp * sizeof(*(phase->dZ_dx)));

    return phase;
}

void flash_calculation_phase_free(PHASE **phase)
{
    PHASE *phase0 = *phase;

    free(phase0->dAx);
    free(phase0->dBx);
    free(phase0->dZ_dx);
    free(phase0->fug);
    free(phase0->phi);
    free(phase0->dphi);
    free(phase0->dphi_dT);
    free(phase0->dphi_dx);
    free(phase0->dfug);
    free(phase0->ai);
    free(phase0->dai_dT);
    free(phase0->dai_dT2);
    free(phase0->da);
    free(phase0->dda_dT);
    free(phase0->dda_dT2);
    free(phase0->bi);
    free(phase0->db);

    free(*phase);
}

/*
# ### Calculate phase parameters
# #### Calculate parameters $a_i$ and $b_i$ of each component
# For PR EOS:
# $$
# b_i = \frac{0.07780 R T_{c,i}}{P_{c,i}}, \\
# a_i = \frac{0.45724 R^2 T^2_{c,i}}{P_{c,i}} (1 + \lambda_i (1 - (\frac{T}{T_{c,i}})^{0.5})) \\
# \lambda_i = 
# \left\{
# \begin{array}{ccc}
# 0.37464 + 1.5432 \omega_i - 0.26992 \omega_i^2 & \text{when} &\omega_i < 0.49\\
# 0.3796 + 1.485 \omega_i - 0.1644 \omega_i^2 + 0.01666 \omega_i^3  & \text{when} & \omega_i >= 0.49
# \end{array}
# \right.
# $$
# 
# For SRK EOS:
# $$
# b_i = \frac{0.08664 R T_{c,i}}{P_{c,i}}, \\
# a_i = \frac{0.42747 R^2 T^2_{c,i}}{P_{c,i}} (1 + \lambda_i (1 - (\frac{T}{T_{c,i}})^{0.5})) \\
# \lambda_i = 0.48 + 1.574 \omega_i - 0.176 \omega_i^2 
# $$
# #### Calculate phase parameters and their derivatives
# $$
# \left\{
# \begin{array}{l}
# a = \sum_{i = 1}^{N_c} x_i S_i \\
# S_i = \sqrt{a_i} \sum_{i = 1}^{N_c} x_j (1 - k_{i,j}) \sqrt{a_j} \\
# b = \sum_{i = 1}^{N_c} x_i b_i
# \end{array}
# \right.
# $$
# and
# $$
# \left\{
# \begin{array}{l}
# \frac{\partial a}{\partial x_i} = 2 S_i, \quad S_i = \sqrt{a_i} \sum_{i = 1}^{N_c} x_j (1 - k_{i,j}) \sqrt{a_j} \\
# \frac{\partial b}{\partial x_i} = b_i
# \end{array}
# \right.
# $$
# where $k_{i,j}$ is the binary interaction parameter.
# $$
# \left\{
# \begin{array}{l}
# A = \frac{a P}{R^2 T^2}\\
# B = \frac{b P}{R T}
# \end{array}
# \right.
# $$
*/

void flash_calculation_compute_phase_parameter(PHASE *phase)
{
    double R, T, P;
    int ncomp, i, j;
    EOS *eos;
    COMP *comp;
    double a, b, da_dT, da_dT2;

    ncomp = phase->ncomp;
    eos = phase->eos;
    comp = eos->comp_list->comp;
    T = eos->temp;
    P = eos->pres;
    R = phase->R;

    if (eos->type == 0) {
        for (i = 0; i < ncomp; i++) {
            double AC, TC, PC, lambda_i, tmp, alpha_i;
            double dtmp, dalpha_i, ddtmp, ddalpha_i;

            AC = eos->comp_list->comp[i].AC;
            TC = eos->comp_list->comp[i].TC;
            PC = eos->comp_list->comp[i].PC;

            lambda_i = 0.48 + 1.574 * AC - 0.176 * AC * AC;

            tmp = 1.0 + lambda_i * (1.0 - sqrt(T / TC));
            alpha_i = tmp * tmp;

            phase->ai[i] = 0.42747 * alpha_i * R * R * TC * TC / PC;
            phase->bi[i] = 0.08664 * R * TC / PC;            

            dtmp = lambda_i * (-1.0) * 0.5 / sqrt(T / TC) / TC;
            dalpha_i = 2.0 * tmp * dtmp;
            phase->dai_dT[i] = 0.42747 * dalpha_i * R * R * TC * TC / PC;

            ddtmp = lambda_i * (-1.0) * 0.5 / TC * (-0.5) * pow(T / TC, -1.5) / TC;
            ddalpha_i = 2.0 * (dtmp * dtmp + tmp * ddtmp);
            phase->dai_dT2[i] = 0.42747 * ddalpha_i * R * R * TC * TC / PC;
        }
    }
    else if (eos->type == 1) {
        for (i = 0; i < ncomp; i++) {
            double AC, TC, PC, lambda_i, tmp, alpha_i;
            double dtmp, dalpha_i, ddtmp, ddalpha_i;

            AC = eos->comp_list->comp[i].AC;
            TC = eos->comp_list->comp[i].TC;
            PC = eos->comp_list->comp[i].PC;

            if (AC < 0.49) {
                lambda_i = 0.37464 + 1.54226 * AC - 0.26992 * AC * AC;
            }
            else {
                lambda_i = 0.3796 + 1.485 * AC - 0.1644 * AC * AC + 0.01666 * AC * AC * AC;
            }

            tmp = 1.0 + lambda_i * (1.0 - sqrt(T / TC));
            alpha_i = 1.0 + lambda_i * (1.0 - sqrt(T / TC));
            alpha_i = alpha_i * alpha_i;

            phase->ai[i] = 0.45724 * alpha_i * R * R * TC * TC / PC;
            phase->bi[i] = 0.077796 * R * TC / PC;

            dtmp = lambda_i * (-1.0) * 0.5 / sqrt(T / TC) / TC;
            dalpha_i = 2.0 * tmp * dtmp;
            phase->dai_dT[i] = 0.45724 * dalpha_i * R * R * TC * TC / PC;

            ddtmp = lambda_i * (-1.0) * 0.5 / TC * (-0.5) * pow(T / TC, -1.5) / TC;
            ddalpha_i = 2.0 * (dtmp * dtmp + tmp * ddtmp);
            phase->dai_dT2[i] = 0.45724 * ddalpha_i * R * R * TC * TC / PC;
        }
    }


    a = 0.0;
    b = 0.0;
    da_dT = 0.0;
    da_dT2 = 0.0;
    for (i = 0; i < ncomp; i++) {
        double da, dda_dT, dda_dT2;

        b += phase->mf[i] * phase->bi[i];
        phase->db[i] = phase->bi[i];

        da = 0.0;
        dda_dT = 0.0;
        dda_dT2 = 0.0;

        for (j = 0; j < ncomp; j++) { 
            a += phase->mf[i] * phase->mf[j] * (1.0 - comp[i].binary[j]) * sqrt(phase->ai[i] * phase->ai[j]);
            da += phase->mf[j] * (1.0 - comp[i].binary[j]) * sqrt(phase->ai[j]);

            da_dT += phase->mf[i] * phase->mf[j] * (1.0 - comp[i].binary[j]) * 0.5 / sqrt(phase->ai[i] * phase->ai[j]) 
                * (phase->ai[j] * phase->dai_dT[i] + phase->ai[i] * phase->dai_dT[j]);
            da_dT2 += phase->mf[i] * phase->mf[j] * (1.0 - comp[i].binary[j]) * 0.5 / sqrt(phase->ai[i] * phase->ai[j]) 
                * (phase->ai[j] * phase->dai_dT2[i] + phase->dai_dT[j] * phase->dai_dT[i] + phase->ai[i] * phase->dai_dT2[j] + phase->dai_dT[i] * phase->dai_dT[j]) 
                + phase->mf[i] * phase->mf[j] * (1.0 - comp[i].binary[j]) * 0.5 * (-0.5) * pow(phase->ai[i] * phase->ai[j], -1.5) 
                * pow(phase->ai[j] * phase->dai_dT[i] + phase->ai[i] * phase->dai_dT[j], 2.0);


            dda_dT += phase->mf[j] * (1.0 - comp[i].binary[j]) * 0.5 / sqrt(phase->ai[j]) * phase->dai_dT[j];
            dda_dT2 += phase->mf[j] * (1.0 - comp[i].binary[j]) * 0.5 * (-0.5) * pow(phase->ai[j], -1.5) * pow(phase->dai_dT[j], 2.0) 
                + phase->mf[j] * (1.0 - comp[i].binary[j]) * 0.5 / sqrt(phase->ai[j]) * phase->dai_dT2[j];
        }

        phase->da[i] = 2.0 * sqrt(phase->ai[i]) * da;

        phase->dda_dT[i] = 2.0 * sqrt(phase->ai[i]) * dda_dT 
            + 2.0 * 0.5 / sqrt(phase->ai[i]) * phase->dai_dT[i] * da;
        phase->dda_dT2[i] = 2.0 * sqrt(phase->ai[i]) * dda_dT2 
            + 2.0 * 0.5 / sqrt(phase->ai[i]) * phase->dai_dT[i] * dda_dT 
            + 2.0 * 0.5 / sqrt(phase->ai[i]) * phase->dai_dT[i] * da_dT 
            + 2.0 * 0.5 / sqrt(phase->ai[i]) * phase->dai_dT2[i] * da 
            + 2.0 * 0.5 * (-0.5) * pow(phase->ai[i], -1.5) * pow(phase->dai_dT[i], 2.0) * da;
    }

    phase->a = a;
    phase->b = b;
    phase->da_dT = da_dT;
    phase->da_dT2 = da_dT2;

    /* A, B */

    phase->A = a * P / (R * R * T * T);
    phase->dAp = phase->A / P;

    phase->dA_dT = da_dT * phase->A / a + a * P / (R * R) * (-2.0) / (T * T * T);
    phase->dA_dT2 = da_dT2 * phase->A / a 
        + da_dT * phase->dA_dT / a
        + da_dT * phase->A * (-1.0) / (a * a) * da_dT 
        + da_dT * P / (R * R) * (-2.0) / (T * T * T)
        + a * P / (R * R) * (-2.0) * (-3.0) / (T * T * T * T);

    for (i = 0; i < ncomp; i++) {
        phase->dAx[i] = phase->A / phase->a * phase->da[i];
    }

    phase->B = b * P / (R * T);
    phase->dBp = phase->B / P;

    phase->dB_dT = b * P / R * (-1.0) / (T * T);
    phase->dB_dT2 = b * P / R * (-1.0) * (-2.0) / (T * T * T);

    for (i = 0; i < ncomp; i++) {
        phase->dBx[i] = phase->B / phase->b * phase->db[i];
    }
}

/*
# ### Calculate compressibility factore Z and its derivatives
# Solve the equation 
# $$Z^3 + c_2 Z^2 + c_1 Z + c_0 = 0,$$
# where 
# $$
# \left\{
# \begin{array}{l}
# c_2 = (u - 1) B - 1 \\
# c_1 = A + (w - u) B^2 - u B \\
# c_0 = - A B - w B^2 - w B^3
# \end{array}
# \right.
# $$
# to get the compressibility factore Z. 
# The derivatives are
# $$
# \left\{
# \begin{array}{l}
# \frac{\partial Z}{\partial P} = - \frac{\frac{\partial c_2}{\partial P} Z^2 
#                                     + \frac{\partial c_1}{\partial P} Z + \frac{\partial c_0}{\partial P}}
#                                 {3 Z^2 + 2 c_2 Z + c_1} \\
# \frac{\partial Z}{\partial x_i} = - \frac{\frac{\partial c_2}{\partial x_i} Z^2 
#                                     + \frac{\partial c_1}{\partial x_i} Z + \frac{\partial c_0}{\partial x_i}}
#                                 {3 Z^2 + 2 c_2 Z + c_1}
# \end{array}
# \right.
# $$
# #### Root selection
# ##### Three real roots
# If three roots are obtained, the middle one is ignored. From the rest two roots, the one that results in the lowest Gibb's free energy will be selected. Let $Z_A$ and $Z_B$ be the two real roots resulting in free energy $G_A$ and $G_B$ respectively. Since free energy 
# $$
# G = \sum_i x_i \log{f_i} \\
# G_A - G_B = \log{\frac{Z_B - B}{Z_A - B}} + \frac{1}{\sigma_2 - \sigma_1} \frac{A}{B} \log{\frac{Z_B + \sigma_2 B}{Z_A + \sigma_2 B} \frac{Z_A + \sigma_1 B}{Z_B + \sigma_1 B}} - (Z_B - Z_A)
# $$
# If $G_A - G_B > 0$, $Z_B$ will be selected and vice versa. For single-phase fluids, if the above scheme selects the largest $Z$ root, the fluid is said to be vapour. Similarly, if the smallest $Z$ root is chosen, the fluid is said to be liquid. 
# ##### One real root
# If only one real root is obtained, we use the method introduced in the paper "**Comparison of Phase Identification Methods Used in Oil Industry Flow Simulations**" to identify phase.
# The thermal expansion coefficent $\alpha$ is defined by:
# $$
# \alpha = (\frac{\partial \log{V}}{\partial T})_P = \frac{1}{V} (\frac{\partial V}{\partial T})_P
# $$
# If $\frac{\partial \alpha}{\partial T} > 0$, the fluid is liquid; otherwise, it is vapour.
*/

void flash_calculation_calculate_compressibility_factor(PHASE *phase)
{
    int ncomp, i, nroot;
    double u, w, A, B, dAp, dBp;
    double c0, c1, c2, dc0, dc1, dc2, 
           *dc0_dx, *dc1_dx, *dc2_dx;
    double dA_dT, dB_dT, dc2_dT, dc1_dT, dc0_dT;
    double dA_dT2, dB_dT2, dc2_dT2, dc1_dT2, dc0_dT2;
    double Q, J, D, Z1, Z2, Z3;
    double Z_h, Z_l;
    double top, down;
    EOS *eos;
    COMP *comp;

    eos = phase->eos;
    comp = eos->comp_list->comp;
    ncomp = phase->ncomp;

    u = eos->para_u;
    w = eos->para_w;
    A = phase->A;
    B = phase->B;
    dAp = phase->dAp;
    dBp = phase->dBp;

    /* Solve the Cubit EOS and get one or three real roots */
    c2 = (u - 1.0) * B - 1.0;
    c1 = A + (w - u) * B * B - u * B;
    c0 = - A * B - w * B * B - w * B * B * B;

    dc2 = (u - 1.0) * dBp;
    dc1 = dAp + (w - u) * 2.0 * B * dBp - u * dBp;
    dc0 = - dAp * B - A * dBp - w * 2.0 * B * dBp - w * 3.0 * B * B * dBp;

    dc2_dx = malloc(ncomp * sizeof(*dc2_dx));
    dc1_dx = malloc(ncomp * sizeof(*dc1_dx));
    dc0_dx = malloc(ncomp * sizeof(*dc0_dx));

    for (i = 0; i < ncomp; i++) {
        dc2_dx[i] = (u - 1.0) * phase->dBx[i];
        dc1_dx[i] = phase->dAx[i] + (w - u) * 2.0 * B * phase->dBx[i] - u * phase->dBx[i];
        dc0_dx[i] = - phase->dAx[i] * B - A * phase->dBx[i] 
            - w * 2.0 * B * phase->dBx[i] - w * 3.0 * B * B * phase->dBx[i];
    }


    dA_dT = phase->dA_dT;
    dB_dT = phase->dB_dT;
    dc2_dT = (u - 1.0) * dB_dT;
    dc1_dT = dA_dT + (w - u) * 2.0 * B * dB_dT - u * dB_dT;
    dc0_dT = - dA_dT * B - A * dB_dT - w * 2.0 * B * dB_dT 
        - w * 3.0 * B * B * dB_dT;

    dA_dT2 = phase->dA_dT2;
    dB_dT2 = phase->dB_dT2;
    dc2_dT2 = (u - 1.0) * dB_dT2;
    dc1_dT2 = dA_dT2 + (w - u) * 2.0 * (dB_dT * dB_dT 
            + B * dB_dT2) - u * dB_dT2;
    dc0_dT2 = - dA_dT2 * B - dA_dT * dB_dT 
        - dA_dT * dB_dT - A * dB_dT2
        - w * 2.0 * (dB_dT * dB_dT + B * dB_dT2) 
        - w * 3.0 * (2.0 * B * dB_dT * dB_dT
                + B * B * dB_dT2);

    Q = (3.0 * c1 - c2 * c2) / 9.0;
    J = (9.0 * c2 * c1 - 27.0 * c0 - 2.0 * c2 * c2 * c2) / 54.0;
    D = Q * Q * Q + J * J;

    Z1 = 0.0;
    Z2 = 0.0;
    Z3 = 0.0;

    if (D > 0.0) {
        double tmp1, tmp2;

        nroot = 1;
        tmp1 = J + sqrt(D);
        tmp2 = J - sqrt(D);

        if (tmp1 > 0) {
            Z1 += pow(tmp1, 1.0 / 3.0);
        }
        else {
            Z1 += - pow(-tmp1, 1.0 / 3.0);
        }

        if (tmp2 > 0) {
            Z1 += pow(tmp2, 1.0 / 3.0);
        }
        else {
            Z1 += - pow(-tmp2, 1.0 / 3.0);
        }

        Z1 += - c2 / 3.0;
    }
    else if (D < 0.0) {
        double theta;

        nroot = 3;
        theta = acos(J / sqrt(- Q * Q * Q));
        Z1 = 2.0 * sqrt(-Q) * cos(theta / 3.0) - c2 / 3.0;
        Z2 = 2.0 * sqrt(-Q) * cos(theta / 3.0 + 2.0 * M_PI / 3.0) - c2 / 3.0;
        Z3 = 2.0 * sqrt(-Q) * cos(theta / 3.0 + 4.0 * M_PI / 3.0) - c2 / 3.0;
    }
    else {
        double tmp;

        nroot = 3;

        if (J > 0) {
            tmp = pow(J, 1.0 / 3.0);
        }
        else {
            tmp = - pow(-J, 1.0 / 3.0);
        }

        Z1 = 2.0 * tmp - c2 / 3.0;
        Z2 = Z3 = - tmp - c2 / 3.0;
    }

    phase->nroot = nroot;

    Z_h = -1e20;
    Z_l = 1e20;
    if (nroot == 3) {
        double sigma_1, sigma_2, dG;

        if (Z_l > Z1) {
            Z_l = Z1;
        }
        if (Z_l > Z2) {
            Z_l = Z2;
        }
        if (Z_l > Z3) {
            Z_l = Z3;
        }

        if (Z_h < Z1) {
            Z_h = Z1;
        }
        if (Z_h < Z2) {
            Z_h = Z2;
        }
        if (Z_h < Z3) {
            Z_h = Z3;
        }

        sigma_1 = eos->para_sigma1;
        sigma_2 = eos->para_sigma2;

        dG = (Z_h - Z_l) + log((Z_l - B) / (Z_h - B))  
            - A / (B * (sigma_2 - sigma_1)) 
            * log((Z_l + sigma_1 * B) / (Z_l + sigma_2 * B) 
                    * (Z_h + sigma_2 * B) / (Z_h + sigma_1 * B));

        if (dG > 0.0) {
            phase->Z = Z_l;
            phase->phase_no = 0;
        }
        else {
            phase->Z = Z_h;
            phase->phase_no = 1;
        }
    }
    else {
        double Tc, Dc;

        Tc = 0;
        Dc = 0;
        for (i = 0; i < ncomp; i++) {
            Tc += phase->mf[i] * comp[i].TC * comp[i].VC;
            Dc += phase->mf[i] * comp[i].VC;
        }
        Tc = Tc / Dc;

        if (eos->temp < Tc) {
            phase->phase_no = 0;
        }
        else {
            phase->phase_no = 1;
        }
        phase->Z = Z1;
    }

    phase->dZ = - (dc2 * phase->Z * phase->Z + dc1 * phase->Z + dc0)         
        / (3.0 * phase->Z * phase->Z + 2.0 * c2 * phase->Z + c1);

    for (i = 0; i < ncomp; i++) {
        phase->dZ_dx[i] = - (dc2_dx[i] * phase->Z * phase->Z 
                + dc1_dx[i] * phase->Z + dc0_dx[i])  
            / (3.0 * phase->Z * phase->Z + 2.0 * c2 * phase->Z + c1);
    }

    phase->dZ_dT = - (dc2_dT * phase->Z * phase->Z + dc1_dT * phase->Z + dc0_dT)    
        / (3.0 * phase->Z * phase->Z + 2.0 * c2 * phase->Z + c1);


    top = - (dc2_dT * phase->Z * phase->Z + dc1_dT * phase->Z + dc0_dT);
    down = 3.0 * phase->Z * phase->Z + 2.0 * c2 * phase->Z + c1;

    phase->dZ_dT2 = - (dc2_dT2 * phase->Z * phase->Z + dc2_dT * 2.0 * phase->Z * phase->dZ_dT   
            + dc1_dT2 * phase->Z + dc1_dT * phase->dZ_dT + dc0_dT2) / down
        + top * (-1.0) / (down * down) * (3.0 * 2.0 * phase->Z * phase->dZ_dT
                + 2.0 * c2 * phase->dZ_dT + 2.0 * dc2_dT * phase->Z
                + dc1_dT);

    if (nroot > 0) {
        double P, Z, T, R, V, dV_dT, dV_dT2, dalpha_dT;

        P = eos->pres;
        T = eos->temp;
        Z = phase->Z;
        R = phase->R;

        V = Z * R * T / P;

        dV_dT = phase->dZ_dT * R * T / P + Z * R / P;
        dV_dT2 = phase->dZ_dT2 * R * T / P + 2.0 * phase->dZ_dT * R / P;

        dalpha_dT = dV_dT2 / V - 1.0 / (V * V) * dV_dT * dV_dT;

        if (dalpha_dT > 0.0) {
            /* liquid */
            phase->phase_no = 0;
        }
        else {
            /* vapor */
            phase->phase_no = 1;
        }
    }

    free(dc2_dx);
    free(dc1_dx);
    free(dc0_dx);
}

/* ### Calculate component fugacity
# $$
# f_i = P x_i \phi_i \\
# \log{\phi_i} = \frac{b_i}{b}(Z - 1) - \log{(Z-B)} 
#             - \frac{A}{B v}(\frac{b_i}{b} - \frac{1}{a} \frac{\partial a}{\partial x_i})
#                 \log{(\frac{2 Z + B(u+v)}{2 Z + B (u-v)})} \\
# v = \sqrt{u^2 - 4 w}
# $$
*/

void flash_calculation_calculate_fugacity(PHASE *phase)
{
    int ncomp = phase->ncomp, i, j;
    double A, B, Z, a, b, u, w, v, p;
    EOS *eos = phase->eos;
    COMP *comp;

    comp = eos->comp_list->comp;

    A = phase->A;
    B = phase->B;
    Z = phase->Z;
    a = phase->a;
    b = phase->b;
    u = eos->para_u;
    w = eos->para_w;
    v = sqrt(u * u - 4.0 * w);
    p = eos->pres;

    for (i = 0; i < ncomp; i++) {
        phase->phi[i] = phase->bi[i] / b * (Z - 1.0) - log(Z - B);
        phase->phi[i] += A / (B * v) * (phase->bi[i] / b - phase->da[i] / a) 
            * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)));
        phase->phi[i] = exp(phase->phi[i]);
        phase->fug[i] = p * phase->mf[i] * phase->phi[i];

        phase->dphi[i] = phase->bi[i] / b * phase->dZ;
        phase->dphi[i] += - (phase->dZ - phase->dBp) / (Z - B);
        phase->dphi[i] += (phase->bi[i] / b - phase->da[i] / a) / v 
            * ((phase->dAp / B - A / (B * B) * phase->dBp)  
                    * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)))  
                    + A / B * ((2.0 * phase->dZ + phase->dBp * (u + v)) / (2.0 * Z + B * (u + v))         
                        - (2.0 * phase->dZ + phase->dBp * (u - v)) / (2.0 * Z + B * (u - v))));
        phase->dphi[i] *= phase->phi[i];
        phase->dfug[i] = phase->fug[i] / p + phase->mf[i] * p * phase->dphi[i];

        /*print("phi[%d]: %lf, pres: %f, mf: %f, phi: %f"%(i, phase.phi[i], phase.cubic_eos->pres, \
          phase.mf[i], phase.phi[i])) */

        for (j = 0; j < ncomp; j++) {
            phase->dphi_dx[i * ncomp + j] = phase->bi[i] / b * phase->dZ_dx[j] 
                - phase->bi[i] / (b * b) * (Z - 1.0) * phase->db[j]
                - (phase->dZ_dx[j]- phase->dBx[j]) / (Z - B)
                + (phase->bi[i] / b - phase->da[i] / a) / v
                * ((phase->dAx[j] / B - phase->dBx[j] * A / (B * B))
                        * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)))
                        + A / B * ((2.0 * phase->dZ_dx[j] + phase->dBx[j] * (u + v))
                            / (2.0 * Z + B * (u + v))
                            - (2.0 * phase->dZ_dx[j] + phase->dBx[j] * (u - v)) 
                            / (2.0 * Z + B * (u - v))))
                + A / (B * v) * log((2.0 * Z + B * (u + v))
                        / (2.0 * Z + B * (u - v)))
                * (- phase->bi[i] / (b * b) * phase->db[j]
                        + phase->da[i] * phase->da[j] / (a * a)
                        - 2.0 / a * sqrt(phase->ai[i] * phase->ai[j]) * (1.0 - comp[i].binary[j]));
            phase->dphi_dx[i * ncomp + j] *= phase->phi[i];
        }

        phase->dphi_dT[i] = phase->bi[i] / b * phase->dZ_dT;
        phase->dphi_dT[i] += - (phase->dZ_dT - phase->dB_dT) / (Z - B);
        phase->dphi_dT[i] += (phase->bi[i] / b - phase->da[i] / a) / v
            * ((phase->dA_dT / B - A / (B * B) * phase->dB_dT) 
                    * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v))) 
                    + A / B * ((2.0 * phase->dZ_dT + phase->dB_dT * (u + v)) / (2.0 * Z + B * (u + v))
                        - (2.0 * phase->dZ_dT + phase->dB_dT * (u - v)) / (2.0 * Z + B * (u - v))));

        phase->dphi_dT[i] += A / B / v * log((2.0 * Z + B * (u + v)) / (2.0 * Z + B * (u - v)))
            * (phase->da[i] / (a * a) * phase->da_dT - phase->dda_dT[i] / a);
        phase->dphi_dT[i] *= phase->phi[i];
    }
}

void flash_calculation_calculate_phase_density(PHASE *phase)
{
    int i, ncomp;
    COMP *comp;
    double mole_den, mole_weight;
    
    mole_den = phase->eos->pres 
        / (phase->Z * phase->R * phase->eos->temp);

    ncomp = phase->ncomp;
    comp = phase->eos->comp_list->comp;
    mole_weight = 0.0;
    for (i = 0; i < ncomp; i++) {
        mole_weight += phase->mf[i] * comp[i].MW;
    }

    phase->density = mole_den * mole_weight;
}





