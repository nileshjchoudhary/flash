#ifndef FC_PHASE
typedef struct PHASE_ {
    int ncomp;
    EOS *eos;
    double *mf;
    
    double R;
    
    double density;
    
    double A;
    double dAp;
    double dA_dT;
    double dA_dT2;
    double *dAx;
    
    double B;
    double dBp;
    double dB_dT;
    double dB_dT2;
    double *dBx;
    
    int nroot;
    double Z;
    double dZ;
    double dZ_dT;
    double dZ_dT2;
    double *dZ_dx;
    
    int phase_no;
    
    double *fug;
    double *phi;
    double *dphi;
    double *dphi_dT;
    double *dphi_dx;
    double *dfug;
    
    double *ai;
    double *dai_dT;
    double *dai_dT2;
    double a;
    double da_dT;
    double da_dT2;
    double *da;
    double *dda_dT;
    double *dda_dT2;
    
    double *bi;
    double b;
    double *db;
} PHASE;

PHASE * flash_calculation_phase_new(EOS *eos, double *mf);
void flash_calculation_phase_free(PHASE **phase);
void flash_calculation_compute_phase_parameter(PHASE *phase);
void flash_calculation_calculate_compressibility_factor(PHASE *phase);
void flash_calculation_calculate_fugacity(PHASE *phase);
void flash_calculation_calculate_phase_density(PHASE *phase);

#define FC_PHASE
#endif

