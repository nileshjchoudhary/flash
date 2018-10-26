#ifndef FC_ANN

typedef enum FLASH_ANN_TRANS_ {
    FLASH_ANN_TRANS_NONE,
    FLASH_ANN_TRANS_LOG,
    FLASH_ANN_TRANS_NEW,

} FLASH_ANN_TRANS;

typedef struct FLASH_TENSOR_ {
    int nr;
    int nc;

    double *value;

} FLASH_TENSOR;

typedef struct FLASH_ANN_ {
    FLASH_TENSOR **W;
    FLASH_TENSOR **b;

    double *feature_scale;
    double *target_scale;
    FLASH_ANN_TRANS trans;

    int level;
} FLASH_ANN;

typedef struct FLASH_STAB_ANN_ {
    FLASH_ANN *ann;

    double safeguard;
    double delta_p;

} FLASH_STAB_ANN;

typedef struct FLASH_SPLIT_ANN_ {
    FLASH_ANN *ann;
    int nK;

} FLASH_SPLIT_ANN;

FLASH_TENSOR * flash_calculation_tensor_read(char *file);
void flash_calculation_tensor_free(FLASH_TENSOR **ft);
FLASH_TENSOR * flash_calculation_tensor_softmax(FLASH_TENSOR *ft, FLASH_TENSOR *result);
FLASH_TENSOR * flash_calculation_tensor_matmul(FLASH_TENSOR *ft1, FLASH_TENSOR *ft2, 
        FLASH_TENSOR *result);
FLASH_TENSOR * flash_calculation_tensor_add(FLASH_TENSOR *ft1, FLASH_TENSOR *ft2,
        FLASH_TENSOR *result);
FLASH_ANN * flash_calculation_ann_model_new(char *file_head, int level, 
        FLASH_ANN_TRANS trans);
void flash_calculation_ann_model_free(FLASH_ANN **ann);
double * flash_calculation_predict_value_with_ANN(FLASH_ANN *ann, 
        FLASH_TENSOR *input);
FLASH_SPLIT_ANN * flash_calculation_split_ann_model_new(char *file_head, int level, 
        char *target_trans, int ncomp);
FLASH_STAB_ANN * flash_calculation_stab_ann_model_new(char *file_head, int level,
        double safeguard, double delta_p);
int flash_calculation_split_ann_predict(FLASH_SPLIT_ANN *fsa, double *input, 
        int n, double *F, double *K);
int flash_calculation_stab_ann_predict(FLASH_STAB_ANN *fsa, double *input,
        int n, int *stable);
void flash_calculation_split_ann_model_free(FLASH_SPLIT_ANN **fsa);
void flash_calculation_stab_ann_model_free(FLASH_STAB_ANN **fsa);

double flash_calculation_split_pred_time_cost(void);
double flash_calculation_stability_pre_time_cost(void);

#define FC_ANN
#endif
