#ifndef FC_DIAGRAM
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

PHASE_DIAGRAM * flash_calculation_phase_diagram_construction(EOS *eos, double *z, 
        double P_start, double P_end, double T_start, double T_end, 
        double *Fv_list, int nF, PHASE_ENVELOPE *pe, CRITICAL_POINT *cp,
        double dP, double dT);
PHASE_DIAGRAM * flash_calculation_draw_phase_diagram(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double *Fv_list, int nF, double dT, double dP, PHASE_ENVELOPE *pe, 
        CRITICAL_POINT *cp, char *output_name);
void flash_calculation_phase_diagram_free(PHASE_DIAGRAM **pd);
#define FC_DIAGRAM
#endif
