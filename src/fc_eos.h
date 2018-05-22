#ifndef FC_EOS
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

EOS * flash_calculation_EOS_new(COMP_LIST *comp_list, double pres, 
        double temp, int type);
#define FC_EOS
#endif
