#ifndef FC_MODEL
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

    int mole;
    int *mole_range;
    int n_mole_range;
    int *mole_dx;
    int n_mole_dx;
    int *mole_d;
    int n_mole_d;

    double dx;

    double *F;
    int nF;

    double T;
    double dxx;
    int selected_component;

    char *ann_trans;

    char *split_ann;
    int split_ann_level;

    char *stab_ann;
    int stab_ann_level;
    double stab_ann_safeguard;
    double stab_ann_delta_p;

    char *output;
    int verb;

} FLASH_MODEL;
FLASH_MODEL * flash_calculation_model_new(char *model_file);
void flash_calculation_flash_model_free(FLASH_MODEL **fm);
#define FC_MODEL
#endif
