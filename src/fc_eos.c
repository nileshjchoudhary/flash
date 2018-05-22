#include "fc.h"

EOS * flash_calculation_EOS_new(COMP_LIST *comp_list, double pres, 
        double temp, int type)
{
    EOS *eos;

    eos = malloc(sizeof(*eos));

    eos->pres = pres;
    eos->temp = temp;
    eos->ncomp = comp_list->ncomp;
    eos->comp_list = comp_list;

    eos->type = type;

    if (type == 0) {
        eos->para_u = 1.0;
        eos->para_w = 0.0;
        eos->para_sigma1 = 1.0;
        eos->para_sigma2 = 0.0;
    }
    else if (type == 1) {
        eos->para_u = 2.0;
        eos->para_w = - 1.0;
        eos->para_sigma1 = 1.0 + sqrt(2.0);
        eos->para_sigma2 = 1.0 - sqrt(2.0);
    }

    return eos;
}

