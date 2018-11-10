#ifndef FC_EOS

typedef enum EOS_TYPE_ {
    EOS_PR = 0,
    EOS_SRK = 1,

} EOS_TYPE;

typedef struct EOS_ {
    double pres;
    double temp;
    
    int ncomp;
    COMP_LIST *comp_list;
    
    EOS_TYPE type;
    double para_u;
    double para_w;
    double para_sigma1;
    double para_sigma2;
} EOS;

EOS * flash_calculation_EOS_new(COMP_LIST *comp_list, double pres, 
        double temp, EOS_TYPE type);
#define FC_EOS
#endif
