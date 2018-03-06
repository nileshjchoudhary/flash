typedef struct CSV_LINE_ {
    int n;
    char **value;
} CSV_LINE;

typedef struct COMP_ {
    char name[20];
    double PC;
    double TC;
    double VC;
    double AC;
    double MW;
    double VS;
    int HYD;
    double *binary;
} COMP;

typedef struct COMP_LIST_ {
    int ncomp;
    COMP *comp;

} COMP_LIST;

typedef struct EOS_ {
    double pres;
    double temp;
    
    int ncomp;
    COMP_LIST *comp_list;
    
    int type;
    double para_u;
    double para_w;
    double para_sigma1;
    double para_sigma2;
} EOS;

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
    
typedef struct CRITICAL_POINT_ {
    double Pc;
    double Tc;
} CRITICAL_POINT;

typedef struct PHASE_LINE_ {
    int n;
    double F;

    double *P;
    double *T;
   
} PHASE_LINE;
    
typedef struct PHASE_DIAGRAM_ {
    PHASE_ENVELOPE *pe;
    CRITICAL_POINT *cp;

    PHASE_LINE **pl;
    int n_line;

} PHASE_DIAGRAM;
    
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
    
    
typedef struct FLASH_TENSOR_ {
    int nr;
    int nc;

    double *value;

} FLASH_TENSOR;

typedef struct FLASH_ANN_ {
    FLASH_TENSOR **W;
    FLASH_TENSOR **b;

    int level;
} FLASH_ANN;

typedef struct FLASH_STAB_ANN_ {
    FLASH_ANN *ann_upper;
    FLASH_ANN *ann_down;

    double safeguard;
    double delta_p;

} FLASH_STAB_ANN;

typedef struct FLASH_SPLIT_ANN_ {
    FLASH_ANN *ann_F;
    FLASH_ANN **ann_K;

    int nK;

} FLASH_SPLIT_ANN;
    
    
typedef struct FLASH_MODEL_ {
    char *type;

    char *prop_file;
    char *bin_file;
    char *z_file;

    double T_min;
    double T_max;
    double P_min;
    double P_max;

    double dP;
    double dT;

    double dx;

    double *F;
    int nF;

    double T;
    double dxx;
    int selected_component;

    char *split_ann;
    int split_ann_level;

    char *stab_ann;
    int stab_ann_level;
    double stab_ann_safeguard;
    double stab_ann_delta_p;

    char *output;
    int verb;

} FLASH_MODEL;
