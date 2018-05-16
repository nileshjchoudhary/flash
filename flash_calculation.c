/* # # Phase Behaviour Calculation Program
# The following code is designed for two phase phase behaviour calculations. The Soave-Redlich-Kwong (SRK) and Peng-Robinson Equation of State are used in this program. The users can choose either of them in their calculations. We provide functions which are capable of doing the following calculations:
# 1. **Phase stability test**
# 2. **Two-phase flash calculation**
# 3. **Saturation pressure or temperature calculations**
# 4. **Critical point calculation**
# 5. **Phase diagram construction**
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include "flash_calculation.h"
#include "mpi.h"
#include <sys/resource.h>
#include <errno.h>
#include <assert.h>

static int verb = 0;
static int split_failure = 0;
static int split_itr = 0;
static double split_solve_time = 0.;
static double split_pred_time = 0.;

static int stab_itr = 0;
static double stab_solve_time = 0.;
static double stab_pred_time = 0.;

static double mem_peak = 0;
static double mem_current = 0;

double flash_calculation_memory_usage(MPI_Comm comm, double *peak)
{
    double a[2], b[2];
    struct rusage RU;

    /* getrusage() */
    getrusage(RUSAGE_SELF, &RU);
    mem_current = RU.ru_maxrss / (double)1024.;

    if (mem_current > mem_peak) mem_peak = mem_current;

    a[0] = mem_current;
    a[1] = mem_peak;
    MPI_Allreduce(a, b, 2, MPI_DOUBLE, MPI_MAX, comm);

    if (peak != NULL) *peak = b[1];

    return b[0];
}

double flash_calculation_get_time(double tarray[])
{
    struct timeval tv;

    if (tarray != NULL) {
        struct rusage RU;

        getrusage(RUSAGE_SELF, &RU);
        tarray[0] = RU.ru_utime.tv_sec + (double)RU.ru_utime.tv_usec * 1e-6;
        tarray[1] = RU.ru_stime.tv_sec + (double)RU.ru_stime.tv_usec * 1e-6;
    }

    /* time */
    gettimeofday(&tv, 0);
    if (tarray != NULL) {
        return tarray[2] = tv.tv_sec + (double)tv.tv_usec * 1e-6;
    }
    else {
        return tv.tv_sec + (double)tv.tv_usec * 1e-6;
    }
}


static CSV_LINE * get_CSV_line(char *line, char *delimiter)
{
    CSV_LINE *csv_line;
    char *tok;
    char **value;
    char delimiter_end[10];

    if (delimiter == NULL) {
        delimiter = " ";
    }
    sprintf(delimiter_end, "%s\n", delimiter);

    csv_line = malloc(sizeof(*csv_line));
    csv_line->n = 0;
    csv_line->value = malloc(10240 * sizeof(*value));

    for (tok = strtok(line, delimiter); tok && *tok; tok = strtok(NULL, delimiter_end))
    {
        int nvalue;

        nvalue = csv_line->n;

        csv_line->value[nvalue] = malloc(10240 * sizeof(char));

        strcpy(csv_line->value[nvalue], tok);
        csv_line->n ++;              
    }

    return csv_line;
}

static void free_CSV_line(CSV_LINE **csv_line)
{
    int i;
    CSV_LINE *csv_line0 = *csv_line;

    for (i = 0; i < csv_line0->n; i++) {
        free(csv_line0->value[i]);
    }

    free(csv_line0->value);
    free(*csv_line);
}

FLASH_TENSOR * flash_calculation_tensor_read(char *file)
{
    FILE *fn;
    char line[10240];
    CSV_LINE **csv_value, *csv_line0;
    int nline, i, j, allocd;
    FLASH_TENSOR *ft;
    char *delimiter = ",";

    allocd = 100;

    csv_value = malloc(100 * sizeof(*csv_value));

    if (!(fn = fopen(file, "r"))) {
        printf("File: %s does not exist!\n", file);
        exit(1);
    }
    nline = 0;
    while (fgets(line, 10240, fn)) {
        csv_line0 = get_CSV_line(line, delimiter);
        csv_value[nline] = csv_line0;

        nline++;

        if (nline >= allocd) {
            csv_value = realloc(csv_value, 
                    (allocd + 100) * sizeof(*csv_value));

            allocd += 100;
        }
    }
    fclose(fn);

    ft = malloc(sizeof(*ft));
    ft->nr = nline;
    ft->nc = csv_value[0]->n;

    ft->value = malloc(ft->nr * ft->nc * sizeof(*(ft->value)));

    for (i = 0; i < ft->nr; i++) {
        for (j = 0; j < ft->nc; j++) {
            ft->value[i * ft->nc + j] = atof(csv_value[i]->value[j]);
        }
    }

    for (i = 0; i < nline; i++) {
        free_CSV_line(&(csv_value[i]));
    }
    free(csv_value);

    return ft;
}

void flash_calculation_tensor_free(FLASH_TENSOR **ft)
{
    FLASH_TENSOR *ft0 = *ft;

    free(ft0->value);
    free(*ft);
}

FLASH_TENSOR * flash_calculation_tensor_softmax(FLASH_TENSOR *ft, FLASH_TENSOR *result)
{
    int i, j;
    double *sum;

    if (result == NULL) {
        result = malloc(sizeof(*result));
        result->nr = ft->nr;
        result->nc = ft->nc;

        result->value = malloc(ft->nr * ft->nc * sizeof(*(result->value)));
    }

    sum = malloc(ft->nr * sizeof(*sum));
    for (i = 0; i < ft->nr; i++) {
        sum[i] = 0.0;
    }

    for (i = 0; i < ft->nr; i++) {
        for (j = 0; j < ft->nc; j++) {
            sum[i] += exp(ft->value[i * ft->nc + j]);
        }
    }

    for (i = 0; i < ft->nr; i++) {
        for (j = 0; j < ft->nc; j++) {
            result->value[i * ft->nc + j] = exp(ft->value[i * ft->nc + j]) / sum[i];
        }
    }

    free(sum);

    return result;
}

FLASH_TENSOR * flash_calculation_tensor_matmul(FLASH_TENSOR *ft1, FLASH_TENSOR *ft2, 
        FLASH_TENSOR *result)
{
    int i, j, k;

    if (ft1->nc != ft2->nr) {
        printf("Tensor 1 dimension: (%d, %d), Tensor 2 dimension: (%d, %d)\n", 
                ft1->nr, ft1->nc, ft2->nr, ft2->nc);
        exit(0);
    }

    if (result == NULL) {
        result = malloc(sizeof(*result));
        result->nr = ft1->nr;
        result->nc = ft2->nc;

        result->value = malloc(result->nr * result->nc * sizeof(*(result->value)));
    }

    for (i = 0; i < result->nr; i++) {
        for (j = 0; j < result->nc; j++) {
            result->value[i * result->nc + j] = 0.0;

            for (k = 0; k < ft1->nc; k++) {
                result->value[i * result->nc + j] 
                    += ft1->value[i * ft1->nc + k] * ft2->value[k * ft2->nc + j];
            }
        }
    }

    return result;
}

FLASH_TENSOR * flash_calculation_tensor_add(FLASH_TENSOR *ft1, FLASH_TENSOR *ft2,
        FLASH_TENSOR *result)
{
    int i, j;

    assert(ft1->nr == ft2->nr);
    assert(ft1->nc == ft2->nc);

    if (result == NULL) {
        result = malloc(sizeof(*result));
        result->nr = ft1->nr;
        result->nc = ft1->nc;

        result->value = malloc(result->nr * result->nc * sizeof(*(result->value)));
    }

    for (i = 0; i < result->nr; i++) {
        for (j = 0; j < result->nc; j++) {
            result->value[i * result->nc + j] = ft1->value[i * ft1->nc + j]
                + ft2->value[i * ft2->nc + j];
        }
    }

    return result;
}

static void flash_calculation_tensor_dump(FLASH_TENSOR *ft)
{
    int i, j;

    printf("Number of rows: %d\n", ft->nr);
    printf("Number of cols: %d\n", ft->nc);
    for (i = 0; i < ft->nr; i++) {
        for (j = 0; j < ft->nc; j++) {
            printf("%lf ", ft->value[i * ft->nc + j]);
        }
        printf("\n");
    }
}

FLASH_ANN * flash_calculation_ann_model_new(char *file_head, int level)
{
    char file[10240];
    int i;
    FLASH_ANN *ann;

    if (level == 0) {
        return NULL;
    }

    ann = malloc(sizeof(*ann));
    ann->level = level;
    ann->W = malloc(level * sizeof(*(ann->W)));
    ann->b = malloc(level * sizeof(*(ann->b)));

    for (i = 0; i < level; i++) {
        sprintf(file, "%s-W%d.csv", file_head, i);

#if 0
        printf("file name: %s\n", file);

#endif
        ann->W[i] = flash_calculation_tensor_read(file);

#if 0
        printf("W[%d]: \n", i);
        flash_calculation_tensor_dump(ann->W[i]);
#endif

        sprintf(file, "%s-b%d.csv", file_head, i);
        ann->b[i] = flash_calculation_tensor_read(file);

#if 0
        printf("b[%d]: \n", i);
        flash_calculation_tensor_dump(ann->b[i]);
#endif
    }

    return ann;
}

void flash_calculation_ann_model_free(FLASH_ANN **ann)
{
    int i;
    FLASH_ANN *ann0 = *ann;

    for (i = 0; i < ann0->level; i++) {
        flash_calculation_tensor_free(&(ann0->W[i]));
        flash_calculation_tensor_free(&(ann0->b[i]));
    }

    free(ann0->W);
    free(ann0->b);

    free(*ann);
}

double flash_calculation_predict_value_with_ANN(FLASH_ANN *ann, 
        FLASH_TENSOR *input)
{
    int level = ann->level, i;
    double value;
    FLASH_TENSOR **y, *x, **W, **b; 

    W = ann->W;
    b = ann->b;

    y = malloc(level * sizeof(*y));

    x = input;
    for (i = 0; i < level - 1; i++) {
        y[i] = flash_calculation_tensor_matmul(x, W[i], NULL);
        flash_calculation_tensor_add(y[i], b[i], y[i]);
        flash_calculation_tensor_softmax(y[i], y[i]);
        x = y[i];
    }

    y[level - 1] = flash_calculation_tensor_matmul(x, W[level - 1], NULL);
    flash_calculation_tensor_add(y[level - 1], b[level - 1], y[level - 1]);
    value = y[level - 1]->value[0];

    for (i = 0; i < level; i++) {
        flash_calculation_tensor_free(&(y[i]));
    }
    free(y);

    return value;
}

FLASH_SPLIT_ANN * flash_calculation_split_ann_model_new(char *file_head, int level, 
        int ncomp)
{
    char file[1024];
    int i;
    FLASH_SPLIT_ANN *fsa;

    printf("    Creating ANN ...\n");
    fsa = malloc(sizeof(*fsa));
    fsa->nK = ncomp;
    fsa->ann_K = malloc(ncomp * sizeof(*(fsa->ann_K)));

    sprintf(file, "%s-Fv", file_head);
    printf("file: %s\n", file);
    fsa->ann_F = flash_calculation_ann_model_new(file, level);

    for (i = 0; i < ncomp; i++) {
        sprintf(file, "%s-K%d", file_head, i);
        printf("file: %s\n", file);

        fsa->ann_K[i] = flash_calculation_ann_model_new(file, level);
    }

    printf("    Done\n");

    return fsa;
}

FLASH_STAB_ANN * flash_calculation_stab_ann_model_new(char *file_head, int level,
        double safeguard, double delta_p)
{
    char file[1024];
    int i;
    FLASH_STAB_ANN *fsa;

    printf("    Creating ANN ...\n");
    fsa = malloc(sizeof(*fsa));

    sprintf(file, "%s-upper", file_head);
    printf("file: %s\n", file);
    fsa->ann_upper = flash_calculation_ann_model_new(file, level);

    sprintf(file, "%s-down", file_head);
    printf("file: %s\n", file);
    fsa->ann_down = flash_calculation_ann_model_new(file, level);

    fsa->safeguard = safeguard;
    fsa->delta_p = delta_p;

    printf("    Done\n");

    return fsa;
}

int flash_calculation_split_ann_predict(FLASH_SPLIT_ANN *fsa, double *input, 
        int n, double *F, double *K)
{
    int i;
    FLASH_TENSOR *ft;

    if (fsa->ann_F == NULL) { return 0; }

    ft = malloc(sizeof(*ft));
    ft->nr = 1;
    ft->nc = n;
    ft->value = malloc(n * sizeof(*(ft->value)));

    for (i = 0; i < n; i++) {
        ft->value[i] = input[i];
    }

    *F = flash_calculation_predict_value_with_ANN(fsa->ann_F, ft);

    for (i = 0; i < fsa->nK; i++) {
        K[i] = flash_calculation_predict_value_with_ANN(fsa->ann_K[i], ft);
    }

    flash_calculation_tensor_free(&ft);

    return 1;
}

int flash_calculation_stab_ann_predict(FLASH_STAB_ANN *fsa, double *input,
        int n, int *stable)
{
    int i;
    FLASH_TENSOR *ft;
    double Pu, Pd, P;

    if (fsa->ann_upper == NULL) {
        return 0;
    }

    ft = malloc(sizeof(*ft));
    ft->nr = 1;
    ft->nc = n;
    ft->value = malloc(n * sizeof(*(ft->value)));

    for (i = 0; i < n - 1; i++) {
        ft->value[i] = input[i];
    }
    P = input[n - 1];

    Pu = flash_calculation_predict_value_with_ANN(fsa->ann_upper, ft);
    Pd = flash_calculation_predict_value_with_ANN(fsa->ann_down, ft);

    flash_calculation_tensor_free(&ft);

    if (((fabs(Pu - P) < fsa->safeguard * fsa->delta_p) 
            || (fabs(Pd - P) < fsa->safeguard * fsa->delta_p))
         && (Pd < Pu)) {
        return 0;
    }

    if (Pu > Pd && P > Pd && P < Pu && Pd > 1.0) {
        *stable = 0;
    }
    else {
        *stable = 1;
    }

    return 1;
}

void flash_calculation_split_ann_model_free(FLASH_SPLIT_ANN **fsa)
{
    int i;
    FLASH_SPLIT_ANN *fsa0 = *fsa;

    if (fsa0 != NULL && fsa0->ann_F !=NULL ) {
        flash_calculation_ann_model_free(&(fsa0->ann_F));
    }

    for (i = 0; i < fsa0->nK; i++) { 
        if (fsa0->ann_K[i] != NULL && fsa0->ann_K[i] != NULL) {
            flash_calculation_ann_model_free(&(fsa0->ann_K[i]));
        }
    }
    free(fsa0->ann_K);

    free(*fsa);
}

void flash_calculation_stab_ann_model_free(FLASH_STAB_ANN **fsa)
{
    int i;
    FLASH_STAB_ANN *fsa0 = *fsa;

    if (fsa0 != NULL && fsa0->ann_upper != NULL) {
        flash_calculation_ann_model_free(&(fsa0->ann_upper));
        flash_calculation_ann_model_free(&(fsa0->ann_down));
    }

    free(fsa0->ann_upper);
    free(fsa0->ann_down);

    free(*fsa);
}

COMP_LIST * flash_calculation_component_new(char *prop_file, char *bin_file)
{
    FILE *fn;
    char line[10240];
    CSV_LINE **csv_value, *csv_line0; 
    int nline, i, j, ncomp;
    COMP_LIST *comp_list;
    COMP *comp;
    double *z;

    csv_value = malloc(20 * sizeof(*csv_value));

    fn = fopen(prop_file, "r");
    nline = 0;
    while (fgets(line, 10240, fn)){
        csv_line0 = get_CSV_line(line, NULL);
        csv_value[nline] = csv_line0;

        nline++;
    }
    fclose(fn);

    ncomp = csv_value[0]->n - 1;
    comp = malloc(ncomp * sizeof(*comp));

    for (i = 0; i < nline; i++) {
        char *name;

        name = csv_value[i]->value[0];

        if (strcmp(name, "PC") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].PC = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "TC") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].TC = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "VC") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].VC = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "AC") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].AC = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "MW") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].MW = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "VS") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].VS = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "HYD") == 0) {
            for (j = 0; j < ncomp; j++) {
                comp[j].HYD = atoi(csv_value[i]->value[j + 1]);
            }
        }
    }

    for (i = 0; i < nline; i++) {
        free_CSV_line(&(csv_value[i]));
    }

    fn = fopen(bin_file, "r");
    nline = 0;
    while (fgets(line, 10240, fn)){
        csv_line0 = get_CSV_line(line, NULL);
        csv_value[nline] = csv_line0;

        nline++;
    }
    fclose(fn);

    for (i = 0; i < nline; i++) {
        comp[i].binary = malloc(ncomp * sizeof(*(comp[i].binary)));

        for (j = 0; j < ncomp; j++) {    
            comp[i].binary[j] = atof(csv_value[i]->value[j]);
        }
    }

    for (i = 0; i < nline; i++) {
        free_CSV_line(&(csv_value[i]));
    }

    comp_list = malloc(sizeof(*comp_list));
    comp_list->ncomp = ncomp;
    comp_list->comp = comp;

    free(csv_value);

    return comp_list;
}

void flash_calculation_component_free(COMP_LIST **comp_list)
{
    int i;
    COMP_LIST *comp_list0 = *comp_list;

    for (i = 0; i < comp_list0->ncomp; i++) {
        free(comp_list0->comp[i].binary);
    }

    free(comp_list0->comp);

    free(*comp_list);
}

double * flash_calculation_composition_new(char *z_file)
{
    FILE *fn;
    char line[10240];
    CSV_LINE **csv_value, *csv_line0; 
    int nline, i, j;
    double *z;

    csv_value = malloc(20 * sizeof(*csv_value));

    fn = fopen(z_file, "r");
    nline = 0;
    while (fgets(line, 10240, fn)){
        csv_line0 = get_CSV_line(line, NULL);
        csv_value[nline] = csv_line0;

        nline++;
    }
    fclose(fn);

    z = malloc(nline * sizeof(*z));

    for (i = 0; i < nline; i++) {
        z[i] = atof(csv_value[i]->value[0]);
    }

    for (i = 0; i < nline; i++) {
        free_CSV_line(&(csv_value[i]));
    }

    free(csv_value);

    return z;
}

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
    int i;
    double mole_den, mole_weight;

    mole_den = phase->eos->pres / (phase->Z * phase->R * phase->eos->temp);

    mole_weight = 0.0;
    for (i = 0; i < phase->ncomp; i++) {
        mole_weight += phase->mf[i] * phase->eos->comp_list->comp[i].MW;
    }

    phase->density = mole_den * mole_weight;
}


/* ## 2. Stability Test
# The following code is designed for stability test. */

/* ### Use Wilson's equation as initial guess for K-values
# $$
# K_i = \frac{P_{c,i}}{P} \exp{5.373 (1 + \omega_i)(1 - \frac{T_{c,i}}{T}})
# $$
*/

double * flash_calculation_estimate_K(EOS *eos, double *K)
{
    int i, ncomp = eos->ncomp;
    double P, T;
    COMP *comp;

    if (K == NULL) {
        K = malloc(ncomp * sizeof(*K));
    }

    P = eos->pres;
    T = eos->temp;
    comp = eos->comp_list->comp;

    for (i = 0; i < ncomp; i++) {
        K[i] = comp[i].PC / P 
            * exp(5.37 * (1.0 + comp[i].AC) * (1.0 - comp[i].TC / T));
    }

    return K;
}

/* ### Composition Initial Guess List
# The initial guess is from the paper "General Strategy for Stability Testing and Phase-Split Calculation in Two and Three Phases" by Zhidong Li and Abbas Firoozabadi in SPE Journal, 2012.
# 
# In total, there are $N_c + 4$ sets of initial guess:
# 1. Wilson's equation: $X_i = K_i x_i$
# 2. Inverse Wilson's equation: $X_i = x_i / K_i$
# 3. $X_i = K_i^{\frac{1}{3}} x_i$
# 4. $X_i = x_i / K_i^{\frac{1}{3}}$
# 5. Pure component: 
# $$
# X_i = 
# \left\{
# \begin{array}[c]
#  0.9 \quad i = j \\
#  0.1 / (N_c - 1) \quad otherwise
# \end{array}
# \right.
# $$ for $j = 1, \cdots, N_c$
*/

double * flash_calculation_stability_analysis_initial_estimate(PHASE *phase)
{
    int ncomp = phase->ncomp, n_guess, i, j;
    double *K, Xi, *est;
    EOS *eos = phase->eos;
    double P, T;

    est = malloc((ncomp + 4) * ncomp * sizeof(*est));
    n_guess = 0;

    P = eos->pres;
    T = eos->temp;

    K = malloc(ncomp * sizeof(*K));
    flash_calculation_estimate_K(eos, K);

    /* Wilson correlation */
    for (i = 0; i < ncomp; i++) {
        Xi = K[i] * phase->mf[i];
        *(est + n_guess * ncomp + i) = Xi;
    }
    n_guess += 1;

    /* Inverse Wilson correlation */
    for (i = 0; i < ncomp; i++) {
        Xi = phase->mf[i] / K[i];
        *(est + n_guess * ncomp + i) = Xi;
    }
    n_guess += 1;

    /* Wilson correlation power(1./3.) */
    for (i = 0; i < ncomp; i++) {
        Xi = pow(K[i], 1.0 / 3.0) * phase->mf[i];
        *(est + n_guess * ncomp + i) = Xi;
    }
    n_guess += 1;

    /* Inverse Wilson correlation power(1./3.) */
    for (i = 0; i < ncomp; i++) {
        Xi = phase->mf[i] / pow(K[i], 1.0 / 3.0); 
        *(est + n_guess * ncomp + i) = Xi;
    }
    n_guess += 1;

    /* A pure phase */
    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                Xi = 0.9;
            }
            else {
                Xi = 0.1 / (ncomp - 1);
            }
            *(est + n_guess * ncomp + j) = Xi;
        }
        n_guess += 1;
    }


    /* A hypothetical idea gas: TODO */
    free(K);

    return est;
}

/* ### Calculate compositions using X
# $$
# x_i = \frac{X_i}{\sum_j X_j} 
# $$
*/

void flash_calculation_calculate_trial_phase_composition(double *X_t, double *x, int ncomp)
{
    double Xs = 0.0;
    int i;

    for (i = 0; i < ncomp; i++) {
        Xs += X_t[i];
    }

    for (i = 0; i < ncomp; i++) {
        x[i] = X_t[i] / Xs;
    }
}

/* ### Update X using SS method
# $$
# X_{t,i} = \exp(\log{x_i} + \log{\phi_i} - \log{\phi_{t,i}})
# $$
*/

void flash_calculation_SS_method_update_X(PHASE *phase, PHASE *phase_t, double *X_t)
{
    int ncomp = phase->ncomp, i;

    for (i = 0; i < ncomp; i++) {
        X_t[i] = exp(log(phase->mf[i]) + log(phase->phi[i]) - log(phase_t->phi[i]));
    }
}

/* ### Calculate derivatives of x with respecte to X
# From $x_i = \frac{X_i}{\sum{X_j}}$, we have
# $$
# \frac{\partial x_i}{\partial X_j} = \frac{\sigma_{i,j}}{\sum{X_j}}
#     - \frac{X_i}{(\sum{X_j})^2}
# $$
# where
# $$
# \sigma_{i,j} = 
# \left\{
# \begin{array}{cc}
# 1 & i = j\\
# 0 & \text{otherwise}
# \end{array}
# \right.
# $$
*/

void flash_calculation_calculate_trial_phase_composition_derivative(double *X_t, double *dx_t, int ncomp)
{
    double sum_X = 0.0, sigma = 0.0;
    int i, j;

    for (i = 0; i < ncomp; i++) {
        sum_X += X_t[i];
    }

    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                sigma = 1.0;
            }
            else {
                sigma = 0.0;
            }

            dx_t[i * ncomp + j] = sigma / sum_X - X_t[i] / (sum_X * sum_X);
        }
    }
}

/* ### Calculate the values of equilibrium equations
# $$
#     D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)}
# $$
*/

void flash_calculation_calculate_stability_equilibrium_equation(PHASE *phase, PHASE *phase_t, 
        double *X_t, double *D)
{
    int ncomp = phase->ncomp, i;

    for (i = 0; i < ncomp; i++) {
        D[i] = log(X_t[i]) - log(phase->mf[i]) + log(phase_t->phi[i]) - log(phase->phi[i]);
    }
}

/* # ### Calculate the derivatives of equilibrium equations with respect to X
# From the equations
# $$
#     D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)},
# $$
# we have the following derivatives
# $$
# \frac{\partial D_i}{\partial X_j} = \frac{\sigma_{i,j}}{X_i} 
#                             + \frac{1}{\phi_i(x)} 
#                             \sum_k{\frac{\partial \phi_i(x)}{\partial x_k} \frac{\partial x_k}{\partial X_j}}
# $$
*/

void flash_calculation_calculate_stability_equilibrium_equation_derivative(PHASE *phase_t, double *dx_t, 
        double *X_t, double *dD)
{
    int ncomp = phase_t->ncomp, i, j, k;
    double sigma = 0.0, temp;

    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                sigma = 1.0;
            }
            else {
                sigma = 0.0;
            }

            dD[i * ncomp + j] = 1.0 / X_t[i] * sigma;

            temp = 0.0;
            for (k = 0; k < ncomp; k++) {
                temp += phase_t->dphi_dx[i * ncomp + k] * dx_t[k * ncomp + j];
            }

            dD[i * ncomp + j] += 1.0 / phase_t->phi[i] * temp;
        }
    }
}

double flash_calculation_calculate_stability_residual(PHASE *phase, PHASE *phase_t, double *X_t, double *res)
{
    int ncomp = phase->ncomp, i;
    double tol = 0.0, tmp1, tmp2;

    for (i = 0; i < ncomp; i++) {
        res[i] = log(X_t[i]) + log(phase_t->phi[i]) - log(phase->mf[i]) - log(phase->phi[i]);

        tmp1 = res[i] * res[i];
        tmp2 = log(phase->mf[i]) + log(phase->phi[i]);
        tmp2 = tmp2 * tmp2;

        if ((tmp1 / tmp2) > tol) {
            tol = tmp1 / tmp2;
        }
    }

    return tol;
}

int flash_calculation_solve_dense_linear_system_LU(int n, double *a0, int pvt[])
{
    int i, j, k;
    double d, r;

#define a(i,j)	(a0[(i) * n + (j)])

    for (i = 0; i < n; i++) {
        d = fabs(a(i,i));
        k = i;
        for (j = i + 1; j < n; j++) {
            if ((r = fabs(a(j,i))) > d) {
                d = r;
                k = j;
            }
        }

        if (d == 0.0) 
            return 0;

        pvt[i] = k;
        if (k != i) { /* exchange row i and row k */
            for (j = i; j < n; j++) {
                d = a(i,j);
                a(i,j) = a(k,j);
                a(k,j) = d;
            }
        }

        if ((d = a(i,i)) != 1.0) {
            a(i,i) = d = 1.0 / d;
            for (j = i + 1; j < n; j++)
                a(i,j) *= d;
        }

        for (j = i + 1; j < n; j++) {
            if ((d = a(j,i)) == 0.0) continue;

            for (k = i + 1; k < n; k++) a(j,k) -= d * a(i,k);
        }
    }

#undef a

    return 1;
}

/* solves L*U*X = B. */
void flash_calculation_solve_dense_linear_system_SV(int n, double *a0, 
        int pvt[], int m, double *b0)
{
    int i, j, k;
    double d;

#define a(i,j)	(a0[(i) * n + (j)])
#define b(i,j)	(b0[(i) * m + (j)])

    for (i = 0; i < n; i++) {
        k = pvt[i];
        if (k != i) {
            /* exchange row i with row k */
            for (j = 0; j < m; j++) {
                d = b(i,j);
                b(i,j) = b(k,j);
                b(k,j) = d;
            }
        }

        if ((d = a(i,i)) != 1.0) {
            for (j = 0; j < m; j++) b(i,j) *= d;
        }

        for (j = i + 1; j < n; j++) {
            if ((d = a(j,i)) == 0.0) continue;
            for (k = 0; k < m; k++) b(j,k) -= d * b(i,k);
        }
    }

    for (i = n - 2; i >= 0; i--)
        for (j = i + 1; j < n; j++) {
            d = a(i,j);
            for (k = 0; k < m; k++) b(i,k) -= d * b(j,k);
        }
#undef a
#undef b
    return;
}

int flash_calculation_solve_dense_linear_system(double *M, double *b, 
        double *x, int n)
{
    int *pvt, i, LU;

    pvt = malloc(n * sizeof(*pvt));

    LU = flash_calculation_solve_dense_linear_system_LU(n, M, pvt);

    if (LU) {
        flash_calculation_solve_dense_linear_system_SV(n, M, 
                pvt, 1, b);

        for (i = 0; i < n; i++) {
            x[i] = b[i];
        }

        free(pvt);
    }
    else {
        free(pvt);

        return LU;
    }

    return 1;
}

/* ### Update X using QNSS method
# Solve $\delta X$ through the following equation
# $$
# J \delta X = - D,
# $$
# where $D_i$ is the value of equilibrium equation, 
# $J$ is Jacobian matrix and $J_{i,j} = \frac{\partial D_i}{\partial X_j}$,
# and update X by 
# $$
# X = X + \delta X
# $$
*/

void flash_calculation_QNSS_method_update_X(double *dD, double *D, double *X_t, int ncomp)
{
    int i;
    double *x;

    x = malloc(ncomp * sizeof(*x));

    for (i = 0; i < ncomp; i++) {
        D[i] = - D[i];
    }

    flash_calculation_solve_dense_linear_system(dD, D, x, ncomp);

    for (i = 0; i < ncomp; i++) {
        X_t[i] += x[i];
    }

    free(x);
}

/* ### Check Stability using X
# 1. If $\sum_i{(\log{\frac{X_i}{z_i}})^2} < \epsilon$, the solution is trivial, try next initial guess.
# 2. If $\sum_i X_i < 1.0$, the phase is stable; otherwise, it is unstable.
*/

int flash_calculation_check_stability(double *X_t, double *z, int ncomp)
{
    int i;
    double sum_Y = 0.0, sum_K = 0.0;

    for (i = 0; i < ncomp; i++) { 
        sum_Y += X_t[i];

        sum_K += log(X_t[i] / z[i]) * log(X_t[i] / z[i]);
    }

    if (sum_K < 1e-2) 
        return -1;

    if (sum_Y <= 1.0 + 1e-8) {
        return 1;
    }
    else {
        return 0;
    }

    return -1;
}

/* ### QNSS method to solve stability analysis
# The tangent-plane criterion of stability analysis of a phase with 
#         compositions z results in solving the following set of equations 
#         (Nghiem and Li, 1984):
# $$
#             D_i = \log{X_i} - \log{z_i} + \log{\phi_i(x)} - \log{\phi_i(z)}
# $$
# where,
# $$
#             x_i = X_i / \sum{X_j}
# $$
# for the primary variables X. 
#         The phase will be stable if (1) $\sum{X_i} < 1$
#                                     (2) $\sum{X_i} = 1$ and $X_i \neq z_i$
#         otherwise the phase is unstable.
#         Equations D are solved using the SS method first and then QNSS 
#         method.
#         To solve the above equations, several sets of initial guesses 
#         should be used to avoid trivial solutions.
*/

int flash_calculation_stability_analysis_QNSS(PHASE *phase, double *K, double tol)
{
    int ncomp = phase->ncomp, i, j, n_guess, itr;
    double *D, *dD, *res, *est;
    double *x_t, *dx_t, *X_t, error;
    PHASE *phase_t;
    int system_status;

    D = malloc(ncomp * sizeof(*D));
    dD = malloc(ncomp * ncomp * sizeof(*dD));
    res = malloc(ncomp * sizeof(*res));

    n_guess = ncomp + 4;
    est = flash_calculation_stability_analysis_initial_estimate(phase);

    x_t = malloc(ncomp * sizeof(*x_t));
    X_t = malloc(ncomp * sizeof(*X_t));
    dx_t = malloc(ncomp * ncomp * sizeof(*dx_t));

    phase_t = flash_calculation_phase_new(phase->eos, x_t);

    flash_calculation_compute_phase_parameter(phase);
    flash_calculation_calculate_compressibility_factor(phase);
    flash_calculation_calculate_fugacity(phase);

    for (i = 0; i < n_guess; i++) {
        for (j = 0; j < ncomp; j++) {
            X_t[j] = est[i * ncomp + j];
        }

        itr = 0;

        while(1) {
            /* Calculate trial phase composition */
            flash_calculation_calculate_trial_phase_composition(X_t, phase_t->mf, ncomp);

            /* Calculation trial phase compostion derivative */
            flash_calculation_calculate_trial_phase_composition_derivative(X_t, dx_t, ncomp);

            /* Calculate the compressibility factor and 
               fugacity coefficient for the trial phase 
               with compositions x_t */
            flash_calculation_compute_phase_parameter(phase_t);
            flash_calculation_calculate_compressibility_factor(phase_t);
            flash_calculation_calculate_fugacity(phase_t);

            /* Calculate the residual */
            error = flash_calculation_calculate_stability_residual(phase, phase_t, X_t, res);

            /* Check if stop */
            if (error < tol)
                break;

            /* Update X_t */
            if (error > 1e-5) {
                flash_calculation_SS_method_update_X(phase, phase_t, X_t);
            }
            else {
                /* Calculate the equilibrim equation and its derivative */
                flash_calculation_calculate_stability_equilibrium_equation(phase, phase_t, X_t, D);
                flash_calculation_calculate_stability_equilibrium_equation_derivative(phase_t, dx_t, X_t, dD);

                /* Update X_t by X_t += - dD^(-1) * D */
                flash_calculation_QNSS_method_update_X(dD, D, X_t, ncomp);
            }

            /* Maximum itrations */
            itr += 1;
            if (itr > 1000) {
                //printf("##### WARNING: Stability_analysis_QNSS reach the maximum iterations!\n");
                break;
            }
        }

        /* Check stability based on Sum(X_t) */
        system_status = flash_calculation_check_stability(X_t, phase->mf, ncomp);

        /* If the solution is trivial, try the next initial guess; 
           otherwise break 
           if system_status is 'Unstable' or system_status is 'Stable':
           break */
        if (system_status != -1) {
            break;
        }
    }

    stab_itr += itr;

    /* if K is not None, we output K = X_t / z */
    if (K != NULL) {
        for (i = 0; i < ncomp; i++) {
            K[i] = X_t[i] / phase->mf[i];
        }
    }

    if (system_status == -1) {
        system_status = 1;
    }


    free(x_t);
    free(X_t);
    free(dx_t);
    free(D);
    free(dD);
    free(res);
    free(est);

    flash_calculation_phase_free(&phase_t);

    return system_status;
}

/* ## 3. Two-phase Flash Calculation
# The following code is designed for two-phase flash calculations.
*/

/* ### Calculate compositions for liquid and vapour phases
# $$
# x_{l,i} = \frac{z_i}{1 + (K_i - 1) F_v} \\
# x_{v,i} = \frac{K_i z_i}{1 + (K_i - 1) F_v}
# $$
# where $z_i$ is the feed composition, $F_v$ is mole fraction of vapour phase.
*/

void flash_calculation_calculate_composition(double *K, double *z, double F_v, 
        double *x_l, double *x_v, int ncomp)
{
    int i;

    for (i = 0; i < ncomp; i++) {
        x_l[i] = z[i] / (1.0 + (K[i] - 1.0) * F_v);
        x_v[i] = K[i] * z[i] / (1.0 + (K[i] - 1.0) * F_v);
    }
}

/* ### Calculate derivatives of compositions for liquid and vapour phases with respect to $K_i$
# $$
# \frac{\partial x_{l,i}}{\partial K_i} = - \frac{z_i F_v}{(1 + (K_i - 1) F_v)^2}  \\
# \frac{\partial x_{v,i}}{\partial K_i} = \frac{z_i}{1 + (K_i - 1) F_v} - \frac{z_i K_i F_v}{(1 + (K_i - 1) F_v)^2}
# $$
*/

void flash_calculation_calculate_composition_derivative(double *K, double *z, double F_v, 
        double *dx_l, double *dx_v, int ncomp)
{
    int i;
    double temp;

    for (i = 0; i < ncomp; i++) {
        temp = 1.0 + (K[i] - 1.0) * F_v;
        dx_l[i] = - z[i] / (temp * temp) * F_v;
        dx_v[i] = z[i] / temp - z[i] * K[i] / (temp * temp) * F_v;
    }
}

/* ### Calculate equilibrium equation
# $$
# G_i = \log{x_{v,i}} + \log{\phi_{v,i}(x_v)} 
#             - \log{x_{l,i}} - \log{\phi_{l,i}(x_l)}
# $$
*/

double flash_calculation_calculate_equilibrium_equation(PHASE *phase_L, PHASE *phase_V, double *G)
{
    int ncomp = phase_L->ncomp, i;
    double error = 0.0;

    for (i = 0; i < ncomp; i++) {
        G[i] = log(phase_V->mf[i]) + log(phase_V->phi[i]);
        G[i] += - log(phase_L->mf[i]) - log(phase_L->phi[i]);

        error += G[i] * G[i];
    }

    return error;
}

/* ### Calculate derivatives of equilibrium equations with respect to K
# $$
# \frac{\partial G_i}{\partial K_j} = \frac{1}{x_{v,i}} \frac{\partial x_{v,i}}{\partial K_i} \sigma_{i,j} 
#                 + \frac{1}{\phi_{v,i}} \frac{\partial \phi_{v,i}}{\partial x_{v,j}} \frac{\partial x_{v,j}}{\partial K_j}
#                - \frac{1}{x_{l,i}} \frac{\partial x_{l,i}}{\partial K_i} \sigma_{i,j} 
#                - \frac{1}{\phi_{l,i}} \frac{\partial \phi_{l,i}}{\partial x_{l,j}} \frac{\partial x_{l,j}}{\partial K_j}     
# $$
*/

void flash_calculation_calculate_equilibrium_equation_derivative(PHASE *phase_L, PHASE *phase_V, 
        double *dx_l, double *dx_v, double *dG)
{
    int ncomp = phase_L->ncomp, i, j;

    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                dG[i * ncomp + j] = 1.0 / phase_V->mf[i] * dx_v[i] - 1.0 / phase_L->mf[i] * dx_l[i];
            }
            else {
                dG[i * ncomp + j] = 0.0;
            }

            dG[i * ncomp + j] += 1.0 / phase_V->phi[i] * phase_V->dphi_dx[i * ncomp + j] * dx_v[j]
                - 1.0 / phase_L->phi[i] * phase_L->dphi_dx[i * ncomp + j] * dx_l[j];
        }
    }
}

/* ### Update K using QNSS method
# Solve
# $$
# dG \delta K = - G,
# $$
# and update $K$ by
# $$
# K = K + \delta K
# $$
*/

void flash_calculation_QNSS_method_update_K(double *dG, double *G, double *K, int ncomp)
{
    int i;
    double *x;

    x = malloc(ncomp * sizeof(*x));

    for (i = 0; i < ncomp; i++) {
        G[i] = -G[i];
    }

    flash_calculation_solve_dense_linear_system(dG, G, x, ncomp);

    for (i = 0; i < ncomp; i++) {
        K[i] += x[i];
    }

    free(x);
}

/* ### Calculate the value of the Rachford-Rice equation
# $$
# V = \sum_i{\frac{z_i (K_i - 1)}{1 + (K_i - 1) F_v}}
# $$
*/

double flash_calculation_calculate_RachfordRice_equation_value(double *K, double *z, double n_V, int ncomp)
{
    int i;
    double value = 0.0;

    for (i = 0; i < ncomp; i++) {
        value += z[i] * (K[i] - 1.0) / (1.0 + (K[i] - 1.0) * n_V);
    }

    return value;
}

/* ### Calculate the derivative of Rachford-Rice equation
# $$
# \frac{\partial V}{\partial F_v} = - \sum_i{\frac{z_i (K_i - 1)^2}{(1 + (K_i - 1)F_v)^2}}
# $$
*/

double flash_calculation_calculate_RachfordRice_equation_derivative(double *K, double *z, double n_V, int ncomp)
{
    int i;
    double value = 0.0;

    for (i = 0; i < ncomp; i++) {
        value += - z[i] * (K[i] - 1.0) * (K[i] - 1.0) 
            / (1.0 + (K[i] - 1.0) * n_V) 
            / (1.0 + (K[i] - 1.0) * n_V);
    }

    return value;
}

/* ### Solve the Rachford-Rice equation
# Solve $F_v$ using an iterative method:
# 1. Solve $\delta F_v$
# $$
# \frac{\partial V}{\partial F_v} \delta F_v = - V 
# $$
# 2. Update $F_v$
# $$
# F_v = F_v + \delta F_v
# $$
# 3. Re-calculate $V$ and $\frac{\partial V}{\partial F_v}$
*/

double flash_calculation_solve_RachfordRice_equation(double *K, double *z, double n_V0, int ncomp)
{
    int i, itr = 0;
    double n_V = n_V0, F, J, d;

    while(1) {
        F = flash_calculation_calculate_RachfordRice_equation_value(K, z, n_V, ncomp);

        if (fabs(F) < 1e-10)
            break;

        J = flash_calculation_calculate_RachfordRice_equation_derivative(K, z, n_V, ncomp);

        d = - F / J;
        n_V += d;

        if (n_V < 0.0) {       
            n_V -= d;
            n_V *= 0.5;
		}

        if (n_V > 1.0) {
            n_V -= d;
            n_V = (1.0 + n_V) * 0.5;
		}
            
        itr += 1;
        if (itr > 100)
            break;
	}       
            
    return n_V;
}

/* ### Update K using SS method
# $$
# K_i = K_i \frac{\phi_{l,i}}{\phi_{v,i}}
# $$
*/
void flash_calculation_SS_method_update_K(double *fug_L, double *fug_V, double *K, int ncomp)
{
	int i;

    for (i = 0; i < ncomp; i++) {
        K[i] = K[i] * fug_L[i] / fug_V[i];
	}
}


/*
    At range [0, 1.0] find a value which makes Rachford-Rich is zero
    using bisection method
*/
double flash_calculation_two_phase_flash_calculation_calculate_initial_F(double *K, double *z, int ncomp) 
{
	double tol = 1e-3, F_min, F_max, 
           V_min, V_max, F_mid, V_mid;
	int itr;
   
    F_min = 0.0;
    F_max = 1.0;
    V_min = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_min, ncomp);
    V_max = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_max, ncomp);
    F_mid = 0.0;
    itr = 0;
    
    while(1) {
        /* calculation the value at bisection position */ 
        F_mid = (F_min + F_max) * 0.5;
        V_mid = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_mid, ncomp);
        
        /* if Value is zero, break */
        if (fabs(V_mid) < 1e-10)
            break;
            
        /* # If V_mid has different signal from V_min,
           #  set F_mid as F_max and V_mid as V_max */
        if (V_mid * V_min < 0.0) {
            V_max = V_mid;
            F_max = F_mid;
		}
        
        /* # If V_mid has different signal from V_max,
           #  set F_mid as F_min and V_mid as V_min */
        if (V_mid * V_max < 0.0) {
            V_min = V_mid;
            F_min = F_mid;
		}
            
        if (fabs(F_min - F_max) < tol) 
            break;
        
        itr += 1;
        if (itr > 100) 
            break;
	}
    
    return F_mid;
}

/* ### QNSS method for Two-phase flash calculation
# Two-phase flahs calculation requires the solution of the following 
# equilibrium and material balance equaions:
# Equilibrium equation:
# $$
#             G_i = \log{x_{v,i}} + \log{\phi_{v,i}(x_v)} 
#             - \log{x_{l,i}} - \log{\phi_{l,i}(x_l)}
# $$
# Material balance equation:
# $$
# V = \sum_i{\frac{z_i (K_i - 1)}{1 + (K_i - 1) F_v}}
# $$
# We take $K_i$ and $F_v$ as the primary variables.
# 
# K-values are solved vis the equilibrium equations, then solve the material balance equation to get the $F_v$. Repeat the process until converge. 
*/
/*
        Two-phase flahs calculation requires the solution of the following 
        equilibrium and material balance equaions:
        Equilibrium equation:
            G[i] = ln(x_v[i]) + ln(phi_v[i]) - ln(x_l[i]) - ln(phi_l[i]) 
            for i = 1,...,Nc
        Material balance equation:
            G[Nc + 1] = Sum_k((z[k] * (K[k] - 1.0)) / (1.0 + F_v * (K[k] - 1.0)))
        We take K[i] and F_v as the primary variables.
        
*/
double flash_calculation_two_phase_flash_Calculation_QNSS(EOS *eos, double *z, 
        double *K, double Fv, double tol)
{
	int i, ncomp = eos->ncomp, itr;
	double F_v, sum_K, error, *x_l, *x_v;
	double *G, *dG, *dx_l, *dx_v, *K0;
    PHASE *phase_L, *phase_V;
    
    /* Initial estimate K */
	K0 = malloc(ncomp * sizeof(*K0));
    if (Fv <= 0.0) {
        flash_calculation_estimate_K(eos, K0);
        F_v = 0.5;

#if 0
        printf("KI: \n");
        for (i = 0; i < ncomp; i++) {
            printf("%e ", K0[i]);
        }
        printf("\n");
#endif
	}
    else {
        sum_K = 0.0;
        for (i = 0; i < ncomp; i++) {
            sum_K += log(K[i]) * log(K[i]);
		}

        if (sum_K < 1e-5) {
            flash_calculation_estimate_K(eos, K0);
            F_v = 0.5;
		}
        else {
			for (i = 0; i < ncomp; i++) 
				K0[i] = K[i];
            F_v = Fv;
		}
	}

#if 0
    {
        double *KK;

        KK = malloc(ncomp * sizeof(*KK));
        flash_calculation_estimate_K(eos, KK);

        printf("Estimate K:\n");
        for (i = 0; i < ncomp; i++) {
            printf("%e ", KK[i]);
        }
        printf("\n");

        free(KK);
    }
#endif

    /* Initial compositions x_v and x_l */
	x_l = malloc(ncomp * sizeof(*x_l));
	x_v = malloc(ncomp * sizeof(*x_v));
    flash_calculation_calculate_composition(K0, z, F_v, x_l, x_v, ncomp);
   
    phase_L = flash_calculation_phase_new(eos, x_l);
    phase_V = flash_calculation_phase_new(eos, x_v);
    
    G = malloc(ncomp * sizeof(*G));
    dG = malloc(ncomp * ncomp * sizeof(*dG));
    dx_l = malloc(ncomp * sizeof(*dx_l));
    dx_v = malloc(ncomp * sizeof(*dx_v));
    
    itr = 0;

    while(1) {
        /* Calculate liquid phase fugacity */
        flash_calculation_compute_phase_parameter(phase_L);
        flash_calculation_calculate_compressibility_factor(phase_L);
        flash_calculation_calculate_fugacity(phase_L);
        
        /* Calculate vapour phase fugacity */
        flash_calculation_compute_phase_parameter(phase_V);
        flash_calculation_calculate_compressibility_factor(phase_V);
        flash_calculation_calculate_fugacity(phase_V);
        
        /* Calculate G error */
        error = flash_calculation_calculate_equilibrium_equation(phase_L, phase_V, G);
        
        /* Check convergence */
        if (error < tol)
            break;
        
        if (error > 1e-5) {
            flash_calculation_SS_method_update_K(phase_L->fug, phase_V->fug, K0, ncomp);
		}
        else {
            /* Calculate the derivatives */
            flash_calculation_calculate_composition_derivative(K0, z, F_v, dx_l, dx_v, ncomp);
            flash_calculation_calculate_equilibrium_equation_derivative(phase_L, phase_V, dx_l, dx_v, dG);
        
            /* Update K */
            flash_calculation_QNSS_method_update_K(dG, G, K0, ncomp);
		}
        
        /* ## Solve Rachford-Rice equation and get F_v */
        F_v = flash_calculation_solve_RachfordRice_equation(K0, z, F_v, ncomp);
        
        /* ## Calculate compositions */
        flash_calculation_calculate_composition(K0, z, F_v, x_l, x_v, ncomp);
        
        itr += 1;
        if (itr > 1000) {
            //printf("##### WARNING: two-phase flash calculation reach maximum iterations!\n");
            break;
		}
	}

    if (verb) {
        printf("Pres: %e, Temp: %e, Itr: %d\n", eos->pres, eos->temp, itr);
    }

    if (fabs(F_v) < 1e-5 || fabs(F_v - 1.0) < 1e-5) {
        split_failure ++;
    }
    split_itr += itr;

    flash_calculation_calculate_phase_density(phase_V);
    flash_calculation_calculate_phase_density(phase_L);

    if (phase_V->density > phase_L->density) {
        F_v = 1.0 - F_v;
        for (i = 0; i < ncomp; i++) {
            K[i] = 1.0 / K0[i];
        }
    }
    else {
        for (i = 0; i < ncomp; i++) {
            K[i] = K0[i];
        }
    }
	
	free(x_l);
	free(x_v);
	free(G);
	free(dG);
	free(dx_l);
	free(dx_v);
	free(K0);

    flash_calculation_phase_free(&phase_L);
    flash_calculation_phase_free(&phase_V);
    /* printf("##### Two-phase flash calculation iterations: %d" %itr); */
    
    return F_v;
}

/* ## 5. Saturation Calculations
# The following code is designed for saturation calculations.
*/

/* ### Search pressure at a given temperature
# At a given temperature, search a pressure $P_s$ at which a single phase exists and 
a pressure $P_u$ at which two phases co-exist. The stability test is used in this function.
*/
void flash_calculation_saturation_calculation_search_Pu_Ps(EOS *eos, double *z, double T, double P_est, 
	int search, double search_step, int Ps_found, int Pu_found, double *Pu_Ps)
{
	/* search: 0 upper
			   1 down */

	int ncomp = eos->ncomp, status;
	double *K, Ps, Pu;
    PHASE *phase;

	eos->temp = T;
    eos->pres = P_est;
    phase = flash_calculation_phase_new(eos, z);
    K = malloc(ncomp * sizeof(*K));
    
    Ps = 0.0;
    Pu = 0.0;
    
    while(1) {
        status = flash_calculation_stability_analysis_QNSS(phase, K, 1e-10);
        
		/* unstable */
        if (status == 0) {
            Pu = eos->pres;
            
            if (search == 0 && !Pu_found) {
                search = 1;
			}
			else if (search == 1 && !Pu_found) {
                search = 0;
			}
                
            Pu_found = 1;
		}
        else {
            Ps = eos->pres;
            Ps_found = 1;
		}
        
        if (Pu_found && Ps_found) 
            break;
        
        if (search == 0)
            eos->pres -= search_step;
        else if (search == 1)
            eos->pres += search_step;
        
        if (eos->pres < 1.0 || eos->pres > 1500.0) {
            Pu = 1.0;
            Ps = 1.0;
            break;
		}
	}

    Pu_Ps[0] = Pu;
	Pu_Ps[1] = Ps;

	free(K);
    flash_calculation_phase_free(&phase);
}

/* ### Search temperature at a given pressure
# At a given pressure, search a temperature $T_s$ at which a single phase exists and a temperature $T_u$ at which two phases co-exist. The stability test is used in this function.
*/
void flash_calculation_saturation_calculation_search_Tu_Ts(EOS *eos, double *z, double P, double T_est, int search, 
	double search_step, int Ts_found, int Tu_found, double *Tu_Ts)
{
	/* default values:  T_est: 100.0, search: 0 upper, 1 down, search_step: 1.0, Ts_found: 0, Tu_found: 0. */
    int ncomp = eos->ncomp;
	double Ts, Tu, *K;
	PHASE *phase;
	int status;

	eos->temp = T_est;
    eos->pres = P;
    phase = flash_calculation_phase_new(eos, z);

    K = malloc(ncomp * sizeof(*K));
    
    Ts = 0.0;
    Tu = 0.0;
    
    while(1) {
        status = flash_calculation_stability_analysis_QNSS(phase, K, 1e-10);
        
        if (status == 0) {
            Tu = eos->temp;
            
            if (search == 0 && !Tu_found) 
                search = 1;
            else if (search == 1 && !Tu_found)
                search = 0;
                
            Tu_found = 1;
		}
        else {
            Ts = eos->temp;
            Ts_found = 1;
		}
        
        if (Tu_found && Ts_found)
            break;
        
        if (search == 0) 
            eos->temp -= search_step;
        else if (search == 1)
            eos->temp += search_step;
        
		if (eos->temp < 1.0 || eos->temp > 1500.0) {
            Tu = 1.0;
            Ts = 1.0;
            break;
		}
	}

	Tu_Ts[0] = Tu;
	Tu_Ts[1] = Ts;

    free(K);
    flash_calculation_phase_free(&phase);
}

/* ### Calculate the pressure at the single-phase and two-phase boundary
# With the $P_s$ and $P_u$, the bisection method is used to search the pressure at the single-phase and two-phase boundary. The method can also be used to generate phase diagram.
*/
double flash_calculation_saturation_calculation_search_boundary_pressure(EOS *eos, double *z, double T, 
	double Pu, double Ps, int search, double tol, int search_status, double *K)
{
	/* search: 0 upper, 1 down */
	/* search_status: -1 anyone, 0 unstable, 1 stable */
	int ncomp = eos->ncomp, status;
	double P_max, P_min, P_mid;
	PHASE *phase;

    if (search == 0) {
        P_max = Ps;
        P_min = Pu;
	}
    else if (search == 1) {
        P_max = Pu;
        P_min = Ps;
	}
    else {
        printf("Saturation search boundary pressure: wrong search!");
	}

    phase = flash_calculation_phase_new(eos, z);
    
    while(1) {
        P_mid = (P_max + P_min) * 0.5;
        
        eos->pres = P_mid;
        
        status = flash_calculation_stability_analysis_QNSS(phase, K, 1e-10);
        
        if (status == 0) {
            if (search == 0) {
                P_min = P_mid;
			}
            else if (search == 1) {
                P_max = P_mid;
			}
		}
        else if (status == 1) {
            if (search == 0)
                P_max = P_mid;
            else if (search == 1)
                P_min = P_mid;
		}
        else {
            printf("Saturation search boundary pressure: wrong stability status!");
		}

        if (fabs(P_max - P_min) < tol) {
            if (search_status == -1 || search_status == status)
                break;
		}
	}

    flash_calculation_phase_free(&phase);

    return eos->pres;
}

/*  ### Calculate Y
# $$
# Y_i = z_i K_i
# $$
*/

void flash_calculation_saturation_calculation_calculate_X(double *K, double *z, double *X, int ncomp)
{
	int i;

    for (i = 0; i < ncomp; i++) {
        X[i] = z[i] * K[i];
	}
}

/* ### Calculate derivative of Y with respect to K
# $$
# \frac{\partial Y_i}{\partial K_i} = z_i
# $$
*/

void flash_calculation_aturation_calculation_calculate_X_derivative(double *z, double *dX_dK, int ncomp)
{
	int i;

    for (i = 0; i < ncomp; i++) {
        dX_dK[i] = z[i];
	}
}

/* ### Calculate fugacity-ratio corrections
# $$
# R_i = \frac{f_{z,i}}{f_{y,i}} (\sum_{j = 1}^{N_c}{Y_i})^{-1}
# $$
*/

void flash_calculation_saturation_calculation_calculate_fugacity_ratio(PHASE *phase, PHASE *phase_t, 
	double *Y, double *Ri)
{
	int i, ncomp = phase->ncomp;
	double sum_Y = 0.0;
    
    for (i = 0; i < ncomp; i++) {
        sum_Y += Y[i];
	}
    
    for (i = 0; i < ncomp; i++) {
        Ri[i] = phase->fug[i] / phase_t->fug[i] / sum_Y;
	}
}

/* ### Update incipient-phase mole numbers with the fungacity-ratio corrections
# $$
# Y_i^{n+1} = Y_i^n R_i
# $$
*/
void flash_calculation_saturation_calculation_update_X_with_fugacity_ratio(double *Y, double *Ri, int ncomp)
{
	int i;
    
    for (i = 0; i < ncomp; i++) {
        Y[i] = Y[i] * Ri[i];
	}
}

/* ### Calculate the Q value
# The recommended approach for determining saturation pressure is based on an approach proposed by Michelsen; he uses the condition:
# $$
# Q(P_{\text{sat}}, y) = 1 - \sum_{i=1}^{N_c}{z_i \frac{\phi_i(z)}{\phi_i(y)}} = 0
# $$
# Then we can have
# $$
# Q(P_{\text{sat}}, y) = 1 - \sum_{i=1}^{N_c}{y_i \frac{f_{z,i}}{f_{y,i}}} = 1 - \sum_{i=1}^{N_c}{Y_i}
# $$
*/
double flash_calculation_saturation_calculation_calculate_Q_value(double *Y, int ncomp)
{
	int i;
    double sum_Y = 0.0;
    
    for (i = 0; i < ncomp; i++) {
        sum_Y += Y[i];
	}
    
    sum_Y = 1.0 - sum_Y;
    
    return sum_Y;
}

/* ### Calculate the derivative of Q with respect to pressure
# $$
# \frac{\partial Q}{\partial P} = \sum_{i=1}^{N_c}{Y_i R_i 
#     (\frac{\partial f_{y,i}}{\partial P} \frac{1}{f_{y,i}}
#     - \frac{\partial f_{z,i}}{\partial P} \frac{1}{f_{z,i}})}
# $$
*/
double flash_calculation_saturation_calculation_calculate_Q_derivative(PHASE *phase, PHASE *phase_t, 
	double *Y, double *Ri)
{
	int i, ncomp = phase->ncomp;
    double dQ = 0.0;
    
    for (i = 0; i < ncomp; i++) {
        dQ += Y[i] * Ri[i] * (phase_t->dfug[i] / phase_t->fug[i] - phase->dfug[i] / phase->fug[i]);
	}

    return dQ;
}

/* ### Saturation Pressure Calculation
# The algorithm is from the book "Phase Behavior" by Curtis H. Whitson and Michael R. Brule.
# 1. Guess a saturation type: bubble- or dewpoint. An incorrect guess will not affect convergence, but the final K values may be "upside down".
# 2. Guess a pressure p*.
# 3. Perform Michelsen's stability test at p*.
# 4. 
#     (a) If the mixture is stable for the current value of p*,this pressure represents p* the upper bound of the search for a saturation 
		pressure on the upper curve of the phase envelope. Return to Step 1 and try a lower pressure to look for an unstable condition.
#     (b) With an unstable condition at p*, this pressure represents the lower bound in the search for a saturation pressure on the upper 
			curve of the phase envelope.
# 5. Having found an unstable solution, use the K values from the stability test to calculate incipient-phase mole numbers at bubble and 
			dewpoint with 
# $$
#                 Y_i = z_i K_i
# $$
# and 
# $$
# Y_i = z_i / K_i
# $$
# 6. Calculate the normalized incipient-phase compositions
# $$
#             y_i = Y_i / \sum{Y_i}
# $$
# 7. Calculate phase Z factors, Z_z and Z_y, and component fugacities, f_z and f_y, from the EOS at the pressure saturation-pressure estimate.
			When multiple Z-factor roots are found for a given phase, the root giving the lowest Gibbs energy should be chosen.
# 8. Calculate fugacity-ratio corrections:
# $$
#               R_i = \frac{f_{z,i}}{f_{y,i}} \sum(Y_i)^(-1)
# $$
# 9. Update incipient-phase mole numbers with the fugacity-ratio corrections:
# $$                    
#                     Y_i = Y_i * R_i^\lambda
# $$
# where four iterations use successive substitution ($\lambda = 1$) followed by a GDEM promotion with lambda given by
# $$
#                     \lambda = \| \frac{b_{11}}{b_{11} - b_{01}}\|
# $$
# where $b_{01} = \sum{\ln{R_i^n} \ln{R_i^{n-1}}}$ and $b_{11} = \sum{ln{R_i^{n-1}} \ln{R_i^{n-1}}}$.
# 10. Calculate a new estimate of saturation pressure using a Newton update:
# $$
#                     P_{\text{sat}}^{n+1} = P_{\text{sat}}^n - \frac{Q^n}{dQ/dp^n}
# $$
# If searching for an upper saturation pressure, the new pressure estimate must be higher than p*. If the new estimate is lower than p*, 
			go to Step 1 and use a new initial-pressure estimate higher than the pressure p* value.
# 11. Check for convergence. Zick suggests the following two criteria (a) 
# $$
# | 1 - \sum{Y_i} | < 10^{-13}
# $$ 
# (b) 
# $$(\sum{\ln{R_i} / \ln{Y_i/z_i}})^2 < 10^{-8}$$
# In addition, check for a trivial solution using the criterion
# $$
#                     \sum{\ln{Y_i/z_i}^2} < 10^{-4}
# $$
# 12. (a) If convergence is not achieved, return to Step 6.
#     (b) If convergence is achieved, determine the saturation type by comparing the mole fraction of the heaviest component in the mixture 
			with that in the incipient phase, where $y_N < z_N$ indicates a bubble point with $K_i = y_i / z_i$ and $y_N > z_N$ indicates a 
			dew point with $K_i = z_i / y_i$, or by comparing the density of the incipient phase with that of the feed.
*/
double flash_calculation_saturation_calculation(EOS *eos, double *z, double T, double P_est, 
        int search, double search_step)
{
	/* search: 0 upper, 1 down */
	int i, itr, ncomp = eos->ncomp;
	double Pu_Ps[2], Pu, Ps, P0;
	double *X, *x, *Ri, *K;
    PHASE *phase_x, *phase_z;

    flash_calculation_saturation_calculation_search_Pu_Ps(eos, z, T, P_est, 
            search, search_step, 0, 0, Pu_Ps);
    Pu = Pu_Ps[0];
	Ps = Pu_Ps[1];

    if (Pu <= 1.0 && Ps <= 1.0) {
        return 1.0;
	}
    
	/* Search unstable pressure */
    K = malloc(ncomp * sizeof(*K));
    P0 = flash_calculation_saturation_calculation_search_boundary_pressure(eos, z, T, 
            Pu, Ps, search, 0.05, 0, K);
    
    if (P0 <= 1.0 && search == 1) {
        return P0;
	}
    
    if (P0 <= 1.0) {
        printf("Cannot find a saturation pressure at a given temperature!");
        return P0;
	}
    
    /* Initial composition */
    X = malloc(ncomp * sizeof(*X));
    x = malloc(ncomp * sizeof(*x));
    Ri = malloc(ncomp * sizeof(*Ri));

    flash_calculation_saturation_calculation_calculate_X(K, z, X, ncomp);
    
    /* Initial phase */
    phase_x = flash_calculation_phase_new(eos, x);
    phase_z = flash_calculation_phase_new(eos, z);
    eos->pres = P0;
    
    itr = 0;
    while(1) { 
		double Q, dQ, dp;

        flash_calculation_calculate_trial_phase_composition(X, x, ncomp);
    
        /* Calculate phase z fugacity */
        flash_calculation_compute_phase_parameter(phase_z);
        flash_calculation_calculate_compressibility_factor(phase_z);
        flash_calculation_calculate_fugacity(phase_z);
        
        /* Calculate phase x fugacity */
        flash_calculation_compute_phase_parameter(phase_x);
        flash_calculation_calculate_compressibility_factor(phase_x);
        flash_calculation_calculate_fugacity(phase_x);
        
        /* Calcualte fugacity-ratio corrections Ri */
        flash_calculation_saturation_calculation_calculate_fugacity_ratio(phase_z, phase_x, 
                X, Ri);
        
		/* # Update incipient-phase mole numbers with
           # fugacity-ratio corrections */
        flash_calculation_saturation_calculation_update_X_with_fugacity_ratio(X, Ri, ncomp);
        
        /* # Calculate a new estimate of saturation pressure
           # using a Newton-Raphson update */
        Q = flash_calculation_saturation_calculation_calculate_Q_value(X, ncomp);
        dQ = flash_calculation_saturation_calculation_calculate_Q_derivative(phase_z, phase_x, X, Ri);
        dp = - Q / dQ;
        eos->pres += dp;

        if (eos->pres < 1.0) {
            eos->pres -= dp;

            eos->pres = (eos->pres + 1.0) / 2.0;
        }
    
        /* Check convergence */
        if (fabs(Q) < 1e-10) 
            break;
            
        itr += 1;
        if (itr > 1000) {
            //printf("##### WARNING: Saturation calculation reach maximum iterations!\n");
            break;
		}
	}
      
    free(K);
    free(X);
    free(x);
    free(Ri);

    flash_calculation_phase_free(&phase_x);
    flash_calculation_phase_free(&phase_z);

    return eos->pres;
}

/* ## 6. Phase Diagram Construction
# The following code is designed for phase diagram construntion.

# ### Calculate the values of equilibrium equations and Rachford-Rice equation
# The equilibrium equations and Rachford-Rice equation are the governing equations
*/
double flash_calculation_phase_diagram_construction_calculate_equations(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *G)
{
    int i, neqn = phase_L->ncomp + 1, ncomp = phase_L->ncomp;
	double error;
    
    /* Equilirium equations */
    flash_calculation_calculate_equilibrium_equation(phase_L, phase_V, G);
    
    /* Material equation */
    G[neqn - 1] = flash_calculation_calculate_RachfordRice_equation_value(K, z, F_v, ncomp);
    
    error = 0.0;
    for (i = 0; i < neqn; i++) {
        error += G[i] * G[i];
	}
    
    return error;
}

/* ### Calculate the derivatives of equilibrium equations and Rachford-Rice equation
# The primary variables are K-values and pressure $P$.
*/
/*
        G[i] = ln(K[i]) + ln(phi_v[i](x_v)) - ln(phi_l[i](x_l))
        dG[i]/dK[j] = 1.0 / K[i] * sigma_ij
                    + 1.0 / phi_v[i](x_v) * dphi_v[i]/dx_v[j] * dx_v[j]/dK[j]
                    - 1.0 / phi_l[i](x_l) * dphi_l[i]/dx_l[i] * dx_l[j]/dK[j]
        dG[i]/dp = 1.0 / phi_v[i](x_v) * dphi_v[i]/dp
                 - 1.0 / phi_l[i](x_l) * dphi_l[i]/dp
        
        G[Nc+1] = Sum_i((K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0)))
        dG[Nc+1]/dK[i] = z[i] / (1.0 + F_v * (K[i] - 1.0))
                    - (K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0))^2 * F_v
        dG[Nc+1]/dp = 0.0
*/
void flash_calculation_phase_diagram_construction_calculate_equations_derivative(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *dG)
{
    int i, j, ncomp = phase_L->ncomp;
    double *dx_l, *dx_v;

    dx_l = malloc(ncomp * sizeof(*dx_l));
    dx_v = malloc(ncomp * sizeof(*dx_v));

    flash_calculation_calculate_composition_derivative(K, z, F_v, dx_l, dx_v, ncomp);

    /* dGi/dKj */
    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                dG[i * (ncomp + 1) + j] = 1.0 / K[i];
			}
            else {
                dG[i * (ncomp + 1) + j] = 0.0;
			}
            
            dG[i * (ncomp + 1) + j] += 1.0 / phase_V->phi[i] * phase_V->dphi_dx[i * ncomp + j] * dx_v[j]
				- 1.0 / phase_L->phi[i] * phase_L->dphi_dx[i * ncomp + j] * dx_l[j];
		}
	}

    /* dGi/dp */
    for (i = 0; i < ncomp; i++) {
		dG[i * (ncomp + 1) + ncomp] = 1.0 / phase_V->phi[i] * phase_V->dphi[i]
			- 1.0 / phase_L->phi[i] * phase_L->dphi[i];
	}
    
    /* dGnc+1/dKi */
    for (i = 0; i < ncomp; i++) {
		double temp;

        temp = 1.0 + F_v * (K[i] - 1.0);
        dG[ncomp * (ncomp + 1) + i] = z[i] / temp 
			- (K[i] - 1.0) * z[i] / (temp * temp) * F_v;
	}
    
    /* dGnc+1/dp */
    dG[ncomp * (ncomp + 1) + ncomp] = 0.0;

	free(dx_l);
	free(dx_v);
}

/*
        G[i] = ln(K[i]) + ln(phi_v[i](x_v)) - ln(phi_l[i](x_l))
        dG[i]/dK[j] = 1.0 / K[i] * sigma_ij
                    + 1.0 / phi_v[i](x_v) * dphi_v[i]/dx_v[j] * dx_v[j]/dK[j]
                    - 1.0 / phi_l[i](x_l) * dphi_l[i]/dx_l[i] * dx_l[j]/dK[j]
        dG[i]/dT = 1.0 / phi_v[i](x_v) * dphi_v[i]/dT
                 - 1.0 / phi_l[i](x_l) * dphi_l[i]/dT
        
        G[Nc+1] = Sum_i((K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0)))
        dG[Nc+1]/dK[i] = z[i] / (1.0 + F_v * (K[i] - 1.0))
                    - (K[i] - 1.0) * z[i] / (1.0 + F_v * (K[i] - 1.0))^2 * F_v
        dG[Nc+1]/dp = 0.0
*/

void flash_calculation_phase_diagram_construction_calculate_equations_derivative_with_K_T(PHASE *phase_L, PHASE *phase_V, 
	double *K, double *z, double F_v, double *dG)
{
    int i, j, ncomp = phase_L->ncomp;
	double *dx_l, *dx_v;
    
    dx_l = malloc(ncomp * sizeof(*dx_l));
    dx_v = malloc(ncomp * sizeof(*dx_v));

    flash_calculation_calculate_composition_derivative(K, z, F_v, dx_l, dx_v, ncomp);

    /* dGi/dKj */
    for (i = 0; i < ncomp; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                dG[i * (ncomp + 1) + j] = 1.0 / K[i];
			}
            else {
                dG[i * (ncomp + 1) + j] = 0.0;
			}
            
            dG[i * (ncomp + 1) + j] += 1.0 / phase_V->phi[i] * phase_V->dphi_dx[i * ncomp + j] * dx_v[j]
				- 1.0 / phase_L->phi[i] * phase_L->dphi_dx[i * ncomp + j] * dx_l[j];
		}
	}
    
    /* dGi/dp */
    for (i = 0; i < ncomp; i++) {
        dG[i * (ncomp + 1) + ncomp] = 1.0 / phase_V->phi[i] * phase_V->dphi_dT[i]
			- 1.0 / phase_L->phi[i] * phase_L->dphi_dT[i];
	}
    
    /* dGnc+1/dKi */
    for (i = 0; i < ncomp; i++) {
		double temp;

        temp = 1.0 + F_v * (K[i] - 1.0);
        dG[ncomp * (ncomp + 1) + i] = z[i] / temp 
			- (K[i] - 1.0) * z[i] / (temp * temp) * F_v;
	}
    
    /* dGnc+1/dp */
    dG[ncomp * (ncomp + 1) + ncomp] = 0.0;

	free(dx_l);
	free(dx_v);
}

double * flash_calculation_phase_diagram_construction_update_variables(double *dG, double *G, int dim)
{
	int i;
	double *x;

	for (i = 0; i < dim; i++) {
		G[i] = -G[i];
	}

	x = malloc(dim * sizeof(*x));

	flash_calculation_solve_dense_linear_system(dG, G, x, dim);
    
    return x;
}

int flash_calculation_search_unstable_temperature(COMP_LIST *comp_list, 
        double *comp_X, double P, double *T_list, int nT)
{
    int i, status;
    EOS *eos;
    PHASE *phase;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);
    phase = flash_calculation_phase_new(eos, comp_X);

    for (i = 0; i < nT; i++) {
        eos->pres = P;
        eos->temp = T_list[i];

        status = flash_calculation_stability_analysis_QNSS(phase, NULL, 1e-10);
        
        if (status == 0) {
            break;
        }
    }

    flash_calculation_phase_free(&phase);
    free(eos);
            
    return i;
}

int flash_calculation_search_stable_temperature(COMP_LIST *comp_list, 
        double *comp_X, double P, double *T_list, int nT)
{
    int i, status;
    EOS *eos;
    PHASE *phase;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);
    phase = flash_calculation_phase_new(eos, comp_X);

    for (i = 0; i < nT; i++) {
        eos->pres = P;
        eos->temp = T_list[i];

        status = flash_calculation_stability_analysis_QNSS(phase, NULL, 1e-10);
        
        if (status == 1) {
            break;
        }
    }
            
    return i;
}


/* ### Saturation Envelope Construction
# We devide the temperature into some intervals wit $\delta T$. At each $T$, we calculate the saturation pressure. 
We first calculate the upper saturation pressure (bubble point pressure), then the lower saturation pressure (dew point pressure).
*/
PHASE_ENVELOPE * flash_calculation_phase_saturation_envelope_construction(EOS *eos, 
        double *z, double T_start, double T_end, double dT, double P_est, double dP)
{
	double dT_min = 0.1, *T_list, P, last;
	int count, n, i, begin_index;
	int flag, flag2, flag3;
    PHASE_ENVELOPE *pe;

	n = (int)((T_end - T_start) / dT) + 1;
	T_list = malloc(n * sizeof(*T_list));

	for (i = 0; i < n - 1; i++) {
		T_list[i] = T_start + dT * i;
	}
	T_list[n - 1] = T_end;

    begin_index = flash_calculation_search_unstable_temperature(eos->comp_list, 
            z, 1.1, T_list, n);

    pe = malloc(sizeof(*pe));
    
    pe->Ps = malloc((2 * n + 20) * sizeof(*(pe->Ps)));
    pe->Ts = malloc((2 * n + 20) * sizeof(*(pe->Ts)));

	count = 0;
    flag = 0;
    P = P_est;

    for (i = 0; i < begin_index; i++) {
        pe->Ps[count] = 1.0;
        pe->Ts[count] = T_list[i];
        count += 1;
    }

    for (i = begin_index; i < n; i++) {
        //printf("----- Upper: temperature %lf\n", T_list[i]);
        P = flash_calculation_saturation_calculation(eos, z, T_list[i], P, 0, dP);
        //printf("             pressure    %lf\n", P); 
        
        if (P > 1.0) {
            pe->Ps[count] = P;
            pe->Ts[count] = T_list[i];
            count += 1;
            flag = 1;
		}
        else {
            if (flag) {
				double T_l, T_r, T_c;

                T_l = T_list[i - 1];
                T_r = T_list[i];
                
                P = pe->Ps[count - 1];
                T_c = (T_l + T_r) * 0.5;

                while(1) {
                    P = flash_calculation_saturation_calculation(eos, z, T_c, P, 0, dP * 0.1);
                    
                    if (P > 1.0) {
                        pe->Ps[count] = P;
                        pe->Ts[count] = T_c;
                        T_l = T_c;
                        count += 1;
					}
                    else {
                        T_r = T_c;
                        P = pe->Ps[count - 1];
					}
                        
                    if (fabs(T_r - T_l) < dT_min) {
                        break;
					}
                    else {
                        T_c = (T_l + T_r) * 0.5;
					}
				}

                break;
			}
            else {
                pe->Ps[count] = 1.0;
                pe->Ts[count] = T_list[i];
                P = P_est;
                count += 1;
			}
		}
	}
    
    flag = 0;
    flag2 = 0;
    flag3 = 0;
    P = pe->Ps[count - 1] - 1.0;
    if (P < 1.0)
        P = 1.0;
	last = count;

    for (i = last - 1; i > 0; i--) {
        //printf("----- Down : temperature %lf\n", pe->Ts[i]);
        if (flag2) {
            pe->Ps[count] = 1.0;
            pe->Ts[count] = pe->Ts[i];
			count += 1;
            continue;
		}

        if (flag3) {
            pe->Ps[count] = pe->Ps[i];
            pe->Ts[count] = pe->Ts[i];
			count += 1;
            continue;
		}
            
        P = flash_calculation_saturation_calculation(eos, z, pe->Ts[i], P, 1, dP * 0.1);
        //printf("             pressure    %lf\n", P);
        
        if (P > 1.0) {
            pe->Ps[count] = P;
            pe->Ts[count] = pe->Ts[i];
            flag = 1;
			count += 1;
		}
        else {
            if (flag) {
                pe->Ps[count] = 1.0;
                pe->Ts[count] = pe->Ts[i];
                flag2 = 1;
				count += 1;
			}
            else {
                P = pe->Ps[i - 1] - 5.0;
                if (P < 1.0)
                    P = 1.0;
            }
		}
	}

    pe->n = count;
    free(T_list);

    return pe;
}

void flash_calculation_phase_envelope_PM_output(PHASE_ENVELOPE_PM *pe_pm,
        int ncomp, int selected_component, char *output)
{
    char file_name[100], file_name_upper[100], file_name_down[100];
    int i, j, k;
    FILE *fp, *fp_upper, *fp_down;
    double xv, P;
    int flag;

    if (output == NULL)
        return;

    sprintf(file_name, "%s-envelope-PM.csv", output);
    fp = fopen(file_name, "a");

    for (k = 0; k < ncomp; k++) {
        fprintf(fp, "Component %d,", k);
    }
    fprintf(fp, "Pressure\n");

    for (i = 0; i < pe_pm->n; i++) {
        for (k = 0; k < ncomp; k++) {
            fprintf(fp, "%lf,", pe_pm->xs[i][k]);
        }

        fprintf(fp, "%lf\n", pe_pm->Ps[i]);
    }
    fclose(fp);

    sprintf(file_name_upper, "%s-envelope-PM-upper.csv", output);
    fp_upper = fopen(file_name_upper, "a");

    sprintf(file_name_down, "%s-envelope-PM-down.csv", output);
    fp_down = fopen(file_name_down, "a");

    flag = 0;
    for (i = 0; i < pe_pm->n; i++) {
        xv = pe_pm->xs[i][selected_component];
        P = pe_pm->Ps[i];

#if 0
        if (i > 1) {
            if (fabs(xv - pe_pm->xs[i - 1][selected_component] > 1e-5)) {
                double dP, dxv, rate, rate0;

                dP = P - pe_pm->Ps[i - 1];
                dxv = xv - pe_pm->xs[i - 1][selected_component];
                rate = dP / dxv;

                dP = pe_pm->Ps[i - 2] - pe_pm->Ps[i - 1];
                dxv = pe_pm->xs[i - 2][selected_component] - pe_pm->xs[i - 1][selected_component];

                if (dxv == 0.0) {
                    rate0 = rate;
                }
                else {
                    rate0 = dP / dxv;
                }

                if (rate > rate0 + 80.0 || rate + 80.0 < rate0) {
                    continue;
                }
            }
        }
#endif

        if (i == 0 || (i > 0 && xv > pe_pm->xs[i - 1][selected_component])
                || (xv == pe_pm->xs[i - 1][selected_component] && flag)) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp_upper, "%lf,", pe_pm->xs[i][j]);
            }

            fprintf(fp_upper, "%lf\n", pe_pm->Ps[i]);

            flag = 0;
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp_down, "%lf,", pe_pm->xs[i][j]);
            }

            fprintf(fp_down, "%lf\n", pe_pm->Ps[i]);

            flag = 1;
        }
    }

    fclose(fp_upper);
    fclose(fp_down);
}

PHASE_ENVELOPE_PM * flash_calculation_phase_saturation_envelope_construction_PM(COMP_LIST *comp_list, 
        double *z, double T, double P_est, double dP, int selected_component, double dx, 
        char *output)
{
    int ncomp = comp_list->ncomp, i, k, n_x_list, count, last;
    double *x_list, *X, sum_no_selected, P;
    PHASE_ENVELOPE_PM *pe_pm;
    int flag, flag2, flag3 = 0;
    EOS *eos;
    FILE *fp;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);

    n_x_list = (int)((1.0 - 0.001) / dx);

    x_list = malloc(n_x_list * sizeof(*x_list));
    for (i = 0; i < n_x_list; i++) {
        x_list[i] = 0.001 + i * dx;
    }

    pe_pm = malloc(sizeof(*pe_pm));
    pe_pm->n = 0;
    pe_pm->Ps = malloc((2 * n_x_list + 20) * sizeof(*(pe_pm->Ps)));
    pe_pm->xs = malloc((2 * n_x_list + 20) * sizeof(*(pe_pm->xs)));

    sum_no_selected = 0.0;
    for (k = 0; k < ncomp; k++) {
        if (k != selected_component) {
            sum_no_selected += z[k];
        }
    }

    count = 0;
    flag = 0;
    P = P_est;

    X = malloc(ncomp * sizeof(*X));

    for (i = 0; i < n_x_list; i++) {
        //printf("X: ");
        for (k = 0; k < ncomp; k++) {
            if (k == selected_component) {
                X[k] = x_list[i];
            }
            else {
                X[k] = z[k] / sum_no_selected * (1.0 - x_list[i]);
            }

            //printf("%lf ", X[k]);
        }
#if 0
        printf("\n");
        printf("T: %lf\n", T);
        printf("Pest: %lf\n", P);
#endif

        P = flash_calculation_saturation_calculation(eos, X, T, P, 0, dP);
        //printf("Result: %lf\n", P);

        if (P > 1.0) {
            pe_pm->Ps[count] = P;
            pe_pm->xs[count] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                pe_pm->xs[count][k] = X[k];
            }
            count += 1;
            flag = 1;
        }
        else {
            if (flag) {
                double x_l, x_r, x_c;

                x_l = x_list[i - 1];
                x_r = x_list[i];

                P = pe_pm->Ps[count - 1];
                x_c = (x_l + x_r) * 0.5;

                for (k = 0; k < ncomp; k++) {
                    if (k == selected_component) {
                        X[k] = x_c;
                    }
                    else {
                        X[k] = z[k] / sum_no_selected * (1.0 - x_c);
                    }
                }

                while(1) {
                    P = flash_calculation_saturation_calculation(eos, X, T, P, 0, dP * 0.1);

                    if (P > 1.0) {
                        pe_pm->Ps[count] = P;
                        pe_pm->xs[count] = malloc(ncomp * sizeof(double));
                        for (k = 0; k < ncomp; k++) {
                            pe_pm->xs[count][k] = X[k];
                        }

                        x_l = x_c;
                        count += 1;
                    }
                    else {
                        x_r = x_c;
                        P = pe_pm->Ps[count - 1];
                    }

                    if (fabs(x_r - x_l) < 0.001) {
                        break;
                    }
                    else {
                        x_c = (x_l + x_r) * 0.5;

                        for (k = 0; k < ncomp; k++) {
                            if (k == selected_component) {
                                X[k] = x_c;
                            }
                            else {
                                X[k] = z[k] / sum_no_selected * (1.0 - x_c);
                            }
                        }
                    }
                }

                break;
            }
            else {
                pe_pm->Ps[count] = 1.0;
                pe_pm->xs[count] = malloc(ncomp * sizeof(double));
                for (k = 0; k < ncomp; k++) {
                    pe_pm->xs[count][k] = X[k];
                }
                count += 1;
            }
        }
    }

    free(X);

    flag = 0;
    flag2 = 0;
    if (pe_pm->Ps[count - 1] == 1.0) {
        flag2 = 1;
        flag3 = 1;
    }

    P = pe_pm->Ps[count - 1] - 1.0;
    if (P < 1.0)
        P = 1.0;
    last = count;

    for (i = last - 1; i >= 0; i--) {
        X = pe_pm->xs[i];

        //printf("----- Down : %d\n", i);
        if (flag2) {
            pe_pm->Ps[count] = 1.0;
            pe_pm->xs[count] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                pe_pm->xs[count][k] = X[k];
            }
            count += 1;
            continue;
        }

        P = flash_calculation_saturation_calculation(eos, X, T, P, 1, dP * 0.1);
        //printf("             pressure    %lf\n", P);

        if (P > 1.0) {
            pe_pm->Ps[count] = P;
            pe_pm->xs[count] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                pe_pm->xs[count][k] = X[k];
            }
            flag = 1;
            count += 1;
        }
        else {
            if (flag) {
                pe_pm->Ps[count] = 1.0;
                pe_pm->xs[count] = malloc(ncomp * sizeof(double));
                for (k = 0; k < ncomp; k++) {
                    pe_pm->xs[count][k] = X[k];
                }
                flag2 = 1;
                count += 1;
			}
            else {
                P = pe_pm->Ps[i - 1] - 5.0;
                if (P < 1.0) {
                    P = 1.0;
                }
            }
		}
    }

    pe_pm->n = count;

    if (flag3) {
        for (i = 0; i < count; i++) {
            free(pe_pm->xs[i]);
        }
        pe_pm->n = 0;
    }

    free(x_list);
    free(eos);

    flash_calculation_phase_envelope_PM_output(pe_pm,
            ncomp, selected_component, output);

    return pe_pm;
}



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
    /* default values:   P_max: 1500.0, P_min: 1.0, P_est: -1.0, K: NULL */
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
    /* default values:   T_max: 1500.0, T_min: 1.0, T_est: -1.0, K: NULL */
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

/* 
# ### Phase Diagram Construction
# The procedure of phase diagram construction is as follows:
# 1. staturation envelope construction
# 2. critical point calculation
# 3. calculate the temperature-pressure line at give $F_v$
*/

PHASE_DIAGRAM * flash_calculation_phase_diagram_construction(EOS *eos, double *z, 
        double P_start, double P_end, double T_start, double T_end, 
        double *Fv_list, int nF, PHASE_ENVELOPE *pe, CRITICAL_POINT *cp, 
        double dP, double dT)
{
    /* default value: dP 1.0, dT 1.0, *pe = NULL */
    /* The number of lines */
    int i, j, ith, ncomp = eos->ncomp, Ts_Ps_count, T_nstep, P_nstep, T_nstep2;
    PHASE_DIAGRAM *pd;
    double Tc, Pc;

    pd = malloc(sizeof(*pd));
    
    /* Calculate critical point */
    if (cp == NULL) {
        pd->cp = flash_calculation_critical_point_calculation(eos, z, 3.5, -1);
    }
    else {
        pd->cp = cp;
    }
    Tc = pd->cp->Tc;
    Pc = pd->cp->Pc;

    /* Construct phase envolope */
    if (pe == NULL) {
        pd->pe = flash_calculation_phase_saturation_envelope_construction(eos, 
                z, T_start, T_end, dT, 250.0, dP);
    }
    else {
        pd->pe = pe;
    }

    pd->pl = malloc(nF * sizeof(*(pd->pl)));
    pd->n_line = nF;
    
    T_nstep = pd->pe->n / 2;

    for (ith = 0; ith < nF; ith++) {
        double Fv, P_last, T_last, P, T;
        double *K, *K0;
        int first;
        double P_start0, T_start0;

        Fv = Fv_list[ith];

        pd->pl[ith] = malloc(sizeof(*(pd->pl[ith])));
        pd->pl[ith]->n = 0;
        pd->pl[ith]->F = Fv;
        pd->pl[ith]->P = malloc(T_nstep * sizeof(*(pd->pl[ith]->P)));
        pd->pl[ith]->T = malloc(T_nstep * sizeof(*(pd->pl[ith]->T)));

        P_last = 0.0;
        T_last = 0.0;
        P = 0.0;
        K = NULL;
        K0 = malloc(ncomp * sizeof(*K0));
        
        /* Increase T and find the corresponding P */
        first = 1;
        for (i = 0; i < T_nstep; i++) { 
            T = pd->pe->Ts[i];

            if (P >= 1.0) {
                K = flash_calculation_calculate_status_at_fixed_temperature_and_Fv(eos, 
                        z, T, Fv, pd->pe->Ps[i], pd->pe->Ps[pd->pe->n - i - 1], P, K);
                P = eos->pres;
            }
            else {
                K = flash_calculation_calculate_status_at_fixed_temperature_and_Fv(eos, 
                        z, T, Fv, pd->pe->Ps[i], pd->pe->Ps[pd->pe->n - i - 1], -1, NULL);
                P = eos->pres;
            }
            
            if (P > 1.0) {
                int k;

                k = pd->pl[ith]->n;

                pd->pl[ith]->P[k] = P;
                pd->pl[ith]->T[k] = T;

                for (j = 0; j < ncomp; j++) {
                    K0[j] = K[j];
                }

                first = 0;

                pd->pl[ith]->n += 1;
            }
            
            if (fabs(P - Pc) < dP && fabs(T - Tc) < dT)
                break;
                
            if (T > Tc && !first)
                break;
        }

        if (K != NULL)  {
            free(K);
            K = NULL;
        }
                
        /* Check the distance between the last found point with critical
        # point. If the pressure distance or the temperature distance 
        # is larger than dP or dT, increase P and find the correspoinding 
        # T
        */
        P_nstep = 0;
        P_start0 = 0.;
        T_start0 = 0.;
        if (pd->pl[ith]->n == 0) {
            P_start0 = P_start;
            T_start0 = T_start;
        }
        else {
            P_start0 = pd->pl[ith]->P[pd->pl[ith]->n - 1];
            T_start0 = pd->pl[ith]->T[pd->pl[ith]->n - 1];
        }
            
        if (fabs(P_start0 - Pc) < dP)
            P_nstep = 0;
        else {
            P_nstep = (int)(fabs(P_start0 - Pc) / dP);
        }

        if (P_nstep > 0) {
            if (P_start0 > Pc)
                dP *= -1.0;
            
            T = 0.0;
            P = P_start0;
            K = NULL;
            
            for (i = 0; i < P_nstep; i++) {
                P += dP;

                if (T > 1.0) {
                    K = flash_calculation_calculate_status_at_fixed_pressure_and_Fv(eos, z, P, Fv,
                            1500.0, 1.0, T, K);
                    T = eos->temp;
                }
                else {
                    K = malloc(ncomp * sizeof(*K));
                    for (j = 0; j < ncomp; j++) {
                        K[j] = K0[j];
                    }

                    K = flash_calculation_calculate_status_at_fixed_pressure_and_Fv(eos, z, P, Fv,
                            1500.0, 1.0, T_start0, K);
                    T = eos->temp;
                }

                if (T > 1.0) {
                    int k;

                    k = pd->pl[ith]->n;

                    pd->pl[ith]->P[k] = P;
                    pd->pl[ith]->T[k] = T;

                    for (j = 0; j < ncomp; j++) {
                        K0[j] = K[j];
                    }

                    P_last = P;
                    T_last = T;

                    pd->pl[ith]->n += 1;
                }
            }

            if (K != NULL) {
                free(K);
                K = NULL;
            }
        }


        T_nstep2 = 0;
        if (fabs(pd->pl[ith]->T[pd->pl[ith]->n - 1] - Tc) < dT) {
            T_nstep2 = 0;
        }
        else {
            T_nstep2 = (int)(fabs(pd->pl[ith]->T[pd->pl[ith]->n - 1] - Tc) / dT);
        }

        if (T_nstep2 > 0) {
            P = 0.0;
            T = pd->pl[ith]->T[pd->pl[ith]->n - 1];
            K = NULL;

            for (i = 0; i < T_nstep2; i++) {
                T -= dT;
                if (P > 1.0) {
                    K = flash_calculation_calculate_status_at_fixed_temperature_and_Fv(eos, z, T, Fv,
                            1500.0, 1.0, P, K);
                    P = eos->pres;
                }
                else {
                    K = malloc(ncomp * sizeof(*K));
                    for (j = 0; j < ncomp; j++) {
                        K[j] = K0[j];
                    }

                    K = flash_calculation_calculate_status_at_fixed_temperature_and_Fv(eos, z, T, Fv, 
                            1500.0, 1.0, (pd->pl[ith]->P[pd->pl[ith]->n - 1] + P_last) * 0.5, K);
                    P = eos->pres;
                }
                
                if (P > 1.0) {
                    int k;

                    k = pd->pl[ith]->n;

                    pd->pl[ith]->P[k] = P;
                    pd->pl[ith]->T[k] = T;

                    pd->pl[ith]->n += 1;
                }
            }

            if (K != NULL) {
                free(K);
                K = NULL;
            }
        }
        
        free(K0);
    }
    
    return pd;
}

/* ## 8. Draw Stability Test Map
# For a feed composition, this function will draw a stability test map which 
    shows the results from stability tests at the given temperature range and 
    pressure range. For the example, the following figure shows the results of 
    stability tests for the temperature in [200 K, 340 K] and the pressure in 
    [20 atm, 80 atm].
*/

void flash_calculation_output_stability_analysis_map(STABILITY_MAP *sm, 
        double *comp_X, int ncomp, char *output_name)
{
    int i, j;
    char file_unstable[100], file_liquid[100], file_vapor[100], file_combination[100];
    FILE *fp;

    sprintf(file_unstable, "%s-unstable.csv", output_name);
    fp = fopen(file_unstable, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure\n");

    for (i = 0; i < sm->n_unstable; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf\n", sm->unstable_temp[i], sm->unstable_pres[i]);
    }
    fclose(fp);

    sprintf(file_liquid, "%s-liquid.csv", output_name);
    fp = fopen(file_liquid, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure\n");

    for (i = 0; i < sm->n_liquid; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf\n", sm->liquid_temp[i], sm->liquid_pres[i]);
    }
    fclose(fp);

    sprintf(file_vapor, "%s-vapor.csv", output_name);
    fp = fopen(file_vapor, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure\n");

    for (i = 0; i < sm->n_vapor; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf\n", sm->vapor_temp[i], sm->vapor_pres[i]);
    }
    fclose(fp);

    sprintf(file_combination, "%s-stability-combination.csv", output_name);
    fp = fopen(file_combination, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure,Unstable,Liquid,Vapor\n");
    for (i = 0; i < sm->n_unstable; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf,1,0,0\n", sm->unstable_temp[i], sm->unstable_pres[i]);
    }

    for (i = 0; i < sm->n_liquid; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf,0,1,0\n", sm->liquid_temp[i], sm->liquid_pres[i]);
    }

    for (i = 0; i < sm->n_vapor; i++) {
        for (j = 0; j < ncomp; j++) {
            fprintf(fp, "%f,", comp_X[j]);
        }
        fprintf(fp, "%lf,%lf,0,0,1\n", sm->vapor_temp[i], sm->vapor_pres[i]);
    }
    fclose(fp);
}

void flash_calculation_output_stability_analysis_map_PM(STABILITY_PM_MAP *sm, 
        double *comp_X, int ncomp, char *output_name)
{
    int i, j;
    char file_unstable[100], file_liquid[100], file_vapor[100], 
         file_combination[100];
    FILE *fp;

    sprintf(file_unstable, "%s-unstable-PM.csv", output_name);
    fp = fopen(file_unstable, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure\n");

    for (i = 0; i < sm->n_unstable; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->unstable_x[i][j]);
            }
        }
        fprintf(fp, "%lf\n", sm->unstable_pres[i]);
    }
    fclose(fp);

    sprintf(file_liquid, "%s-liquid-PM.csv", output_name);
    fp = fopen(file_liquid, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure\n");

    for (i = 0; i < sm->n_liquid; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->liquid_x[i][j]);
            }
        }
        fprintf(fp, "%lf\n", sm->liquid_pres[i]);
    }
    fclose(fp);

    sprintf(file_vapor, "%s-vapor-PM.csv", output_name);
    fp = fopen(file_vapor, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure\n");

    for (i = 0; i < sm->n_vapor; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->vapor_x[i][j]);
            }
        }
        fprintf(fp, "%lf\n", sm->vapor_pres[i]);
    }
    fclose(fp);

    sprintf(file_combination, "%s-stability-combination-PM.csv", output_name);
    fp = fopen(file_combination, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure,Unstable,Liquid,Vapor\n");
    for (i = 0; i < sm->n_unstable; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->unstable_x[i][j]);
            }
        }
        fprintf(fp, "%lf,1,0,0\n", sm->unstable_pres[i]);
    }

    for (i = 0; i < sm->n_liquid; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->liquid_x[i][j]);
            }
        }
        fprintf(fp, "%lf,0,1,0\n", sm->liquid_pres[i]);
    }

    for (i = 0; i < sm->n_vapor; i++) {
        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%f,", sm->vapor_x[i][j]);
            }
        }
        fprintf(fp, "%lf,0,0,1\n", sm->vapor_pres[i]);
    }
    fclose(fp);
}

STABILITY_MAP * flash_calculation_draw_stability_analysis_map(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double dT, double dP, FLASH_STAB_ANN *fsa, char *output_name)
{
    int n_pres, n_temp, i, j, ncomp = comp_list->ncomp;
    double *pres_list, *temp_list;
    EOS *eos;
    PHASE *phase;
    STABILITY_MAP *sm;
    int status;
    char file_unstable[100], file_liquid[100], file_vapor[100], file_combination[100];
    FILE *fp;
    double solve_time, pred_time;

    n_pres = (int)((P_max - P_min) / dP);
    n_temp = (int)((T_max - T_min) / dT);

    pres_list = malloc(n_pres * sizeof(*pres_list));
    temp_list = malloc(n_temp * sizeof(*temp_list));

    for (i = 0; i < n_pres; i++) {
        pres_list[i] = P_min + i * dP;
    }

    for (i = 0; i < n_temp; i++) {
        temp_list[i] = T_min + i * dT;
    }

    sm = malloc(sizeof(*sm));
    sm->n_unstable = 0;
    sm->unstable_pres = malloc(n_pres * n_temp * sizeof(*(sm->unstable_pres)));
    sm->unstable_temp = malloc(n_pres * n_temp * sizeof(*(sm->unstable_temp)));
    sm->n_liquid = 0;
    sm->liquid_pres = malloc(n_pres * n_temp * sizeof(*(sm->liquid_pres)));
    sm->liquid_temp = malloc(n_pres * n_temp * sizeof(*(sm->liquid_temp)));
    sm->n_vapor = 0;
    sm->vapor_pres = malloc(n_pres * n_temp * sizeof(*(sm->vapor_pres)));
    sm->vapor_temp = malloc(n_pres * n_temp * sizeof(*(sm->vapor_temp)));

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);
    phase = flash_calculation_phase_new(eos, comp_X);

    for (i = 0; i < n_pres; i++) {
        for (j = 0; j < n_temp; j++) {
            int flag = 0;

            eos->pres = pres_list[i];
            eos->temp = temp_list[j];

            solve_time = flash_calculation_get_time(NULL);
            if (fsa != NULL) {
                double *input;
                int n, k;

                n = ncomp + 2;
                input = malloc(n * sizeof(*input));

                for (k = 0; k < ncomp; k++) {
                    input[k] = comp_X[k];
                }
                input[ncomp] = temp_list[j];
                input[ncomp + 1] = pres_list[j];

                flag = flash_calculation_stab_ann_predict(fsa, input, n, &status);
            }

            if (!flag) {
                status = flash_calculation_stability_analysis_QNSS(phase, NULL, 1e-10);
            }

            solve_time = flash_calculation_get_time(NULL) - solve_time;
            stab_solve_time += solve_time;

            if (status == 1) {
                if (phase->phase_no == 0) {
                    int k;

                    k = sm->n_liquid;
                    sm->liquid_temp[k] = temp_list[j];
                    sm->liquid_pres[k] = pres_list[i];

                    sm->n_liquid += 1;
                }
                else {
                    int k;

                    k = sm->n_vapor;
                    sm->vapor_temp[k] = temp_list[j];
                    sm->vapor_pres[k] = pres_list[i];

                    sm->n_vapor += 1;
                }
            }
            else if (status == 0) {
                int k;

                k = sm->n_unstable;
                sm->unstable_temp[k] = temp_list[j];
                sm->unstable_pres[k] = pres_list[i];

                sm->n_unstable += 1;
            }
        }
    }

    if (output_name != NULL) {
        flash_calculation_output_stability_analysis_map(sm, 
                comp_X, ncomp, output_name);
    }

    free(pres_list);
    free(temp_list);
    flash_calculation_phase_free(&phase);

    free(eos);

    return sm;
}

STABILITY_PM_MAP * flash_calculation_draw_stability_analysis_map_PM(COMP_LIST *comp_list, 
        double *comp_X, double T, double P_min, double P_max, double dP, int selected_component, 
        double dxx, FLASH_STAB_ANN *fsa, char *output_name)
{
    int n_pres, n_x, i, j, k, ncomp = comp_list->ncomp;
    double *pres_list, *x_list, *x, sum_no_selected;
    EOS *eos;
    PHASE *phase;
    STABILITY_PM_MAP *sm;
    int status;
    char file_unstable[100], file_liquid[100], file_vapor[100], file_combination[100];
    FILE *fp;
    double solve_time, pred_time;

    x = malloc(ncomp * sizeof(*x));
    n_pres = (int)((P_max - P_min) / dP);
    pres_list = malloc(n_pres * sizeof(*pres_list));
    for (i = 0; i < n_pres; i++) {
        pres_list[i] = P_min + i * dP;
    }

    n_x = (int)((1.0 - 0.001) / dxx) + 1;
    x_list = malloc(n_x * sizeof(*x_list));
    for (i = 0; i < n_x - 1; i++) {
        x_list[i] = 0.001 + i * dxx;
    }
    x_list[n_x - 1] = 0.999;

    sm = malloc(sizeof(*sm));
    sm->n_unstable = 0;
    sm->unstable_pres = malloc(n_pres * n_x * sizeof(*(sm->unstable_pres)));
    sm->unstable_x = malloc(n_pres * n_x * sizeof(*(sm->unstable_x)));
    sm->n_liquid = 0;
    sm->liquid_pres = malloc(n_pres * n_x * sizeof(*(sm->liquid_pres)));
    sm->liquid_x = malloc(n_pres * n_x * sizeof(*(sm->liquid_x)));
    sm->n_vapor = 0;
    sm->vapor_pres = malloc(n_pres * n_x * sizeof(*(sm->vapor_pres)));
    sm->vapor_x = malloc(n_pres * n_x * sizeof(*(sm->vapor_x)));

    for (i = 0; i < n_pres * n_x; i++) {
        sm->unstable_x[i] = malloc(ncomp * sizeof(double));
        sm->liquid_x[i] = malloc(ncomp * sizeof(double));
        sm->vapor_x[i] = malloc(ncomp * sizeof(double));
    }

    eos = flash_calculation_EOS_new(comp_list, 0.0, T, 0);
    phase = flash_calculation_phase_new(eos, x);

    sum_no_selected = 0.0;
    for (i = 0; i < ncomp; i++) {
        if (i != selected_component) {
            sum_no_selected += comp_X[i];
        }
    }

    for (i = 0; i < n_pres; i++) {
        for (j = 0; j < n_x; j++) {
            int flag = 0;

            for (k = 0; k < ncomp; k++) {
                if (k == selected_component) {
                    x[k] = x_list[j];
                }
                else {
                    x[k] = (1.0 - x_list[j]) * comp_X[k] / sum_no_selected;
                }
            }

            eos->pres = pres_list[i];
            eos->temp = T;

            solve_time = flash_calculation_get_time(NULL);

            if (fsa != NULL) {
                double *input;
                int n, k;

                n = ncomp + 1;
                input = malloc(n * sizeof(*input));

                for (k = 0; k < ncomp; k++) {
                    input[k] = x[k];
                }
                input[ncomp] = pres_list[i];

                flag = flash_calculation_stab_ann_predict(fsa, input, n, &status);
            }

            if (!flag) {
                status = flash_calculation_stability_analysis_QNSS(phase, NULL, 1e-10);
            }

            solve_time = flash_calculation_get_time(NULL) - solve_time;
            stab_solve_time += solve_time;

            if (status == 1) {
                if (phase->phase_no == 0) {
                    int l;

                    l = sm->n_liquid;
                    sm->liquid_pres[l] = pres_list[i]; 

                    for (k = 0; k < ncomp; k++) {
                        sm->liquid_x[l][k] = x[k];
                    }

                    sm->n_liquid += 1;
                }
                else {
                    int l;

                    l = sm->n_vapor;
                    sm->vapor_pres[l] = pres_list[i];

                    for (k = 0; k < ncomp; k++) {
                        sm->vapor_x[l][k] = x[k];
                    }

                    sm->n_vapor += 1;
                }
            }
            else if (status == 0) {
                int l;

                l = sm->n_unstable;
                sm->unstable_pres[l] = pres_list[i];

                for (k = 0; k < ncomp; k++) {
                    sm->unstable_x[l][k] = x[k];
                }

                sm->n_unstable += 1;
            }
        }
    }

    if (output_name != NULL) {
        flash_calculation_output_stability_analysis_map_PM(sm, 
                NULL, ncomp, output_name);
    }

    free(pres_list);
    free(x_list);
    free(x);
    flash_calculation_phase_free(&phase);

    free(eos);

    return sm;
}

/* ## 8. Draw Two-phase Flash Calculation Map
# For a feed composition, this function will draw a two-phase flash calculation 
map which shows mole fraction of vapour phase at the given temperature range 
and pressure range. For the example, the following figure shows the results 
of two-phase flash calculation for the temperature in [200 K, 340 K] and the 
pressure in [20 atm, 80 atm].
*/

void flash_calculation_output_split_calculation_map(SPLIT_MAP *sm, 
        double *comp_X, int ncomp, int filter, char *output_name)
{
    int i, j;
    char file_name[100];
    FILE *fp;

    sprintf(file_name, "%s-split-calculation.csv", output_name);
    fp = fopen(file_name, "a");
    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure,Fv");
    for (j = 0; j < ncomp; j++) {
        fprintf(fp, ",K_%d", j);
    }
    fprintf(fp, "\n");

    for (i = 0; i < sm->n; i++) {
        if (filter && (fabs(sm->F[i]) < 1e-5 || fabs(sm->F[i] - 1.0) < 1e-5)) {
            continue;
        }

        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", sm->x[i][j]);
            }
        }

        fprintf(fp, "%lf,%lf,%lf", sm->temp[i], sm->pres[i], sm->F[i]);

        for (j = 0; j < ncomp; j++) {
            fprintf(fp, ",%lf", sm->K[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

void flash_calculation_output_split_calculation_map_PM(SPLIT_PM_MAP *sm, 
        double *comp_X, int ncomp, int filter, char *output_name)
{
    int i, j;
    char file_name[100];
    FILE *fp;

    sprintf(file_name, "%s-split-calculation-PM.csv", output_name);
    fp = fopen(file_name, "a");
    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Pressure,Fv");
    for (j = 0; j < ncomp; j++) {
        fprintf(fp, ",K_%d", j);
    }
    fprintf(fp, "\n");

    for (i = 0; i < sm->n; i++) {
        if (filter && (fabs(sm->F[i]) < 1e-5 || fabs(sm->F[i] - 1.0) < 1e-5)) {
            continue;
        }

        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", sm->x[i][j]);
            }
        }

        fprintf(fp, "%lf,%lf", sm->pres[i], sm->F[i]);

        for (j = 0; j < ncomp; j++) {
            fprintf(fp, ",%lf", sm->K[i][j]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
}

SPLIT_MAP * flash_calculation_draw_split_calculation_map(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double dT, double dP, FLASH_SPLIT_ANN *fsa, char *output_name)
{
    int i, j, k, ncomp = comp_list->ncomp;
    int n_pres, n_temp;
    double *pres_list, *temp_list;
    SPLIT_MAP *sm;
    int first, status;
    double Fv, Fv0, *K0, *K, solve_time;
    EOS *eos;
    PHASE *phase;
    char file_name[100];
    FILE *fp;

    n_pres = (int)((P_max - P_min) / dP);
    n_temp = (int)((T_max - T_min) / dT);

    pres_list = malloc(n_pres * sizeof(*pres_list));
    temp_list = malloc(n_temp * sizeof(*temp_list));

    for (i = 0; i < n_pres; i++) {
        pres_list[i] = P_min + i * dP;
    }

    for (i = 0; i < n_temp; i++) {
        temp_list[i] = T_min + i * dT;
    }

    sm = malloc(sizeof(*sm));
    sm->n = n_pres * n_temp;
    sm->temp = malloc(n_pres * n_temp * sizeof(*(sm->temp)));
    sm->pres = malloc(n_pres * n_temp * sizeof(*(sm->pres)));
    sm->F = malloc(n_pres * n_temp * sizeof(*(sm->F)));
    sm->K = malloc(n_pres * n_temp * sizeof(*(sm->K)));
    sm->x = malloc(n_pres * n_temp * sizeof(*(sm->x)));
    for (i = 0; i < n_pres * n_temp; i++) {
        sm->K[i] = malloc(ncomp * sizeof(*(sm->K[i])));
        sm->x[i] = malloc(ncomp * sizeof(*(sm->x[i])));
    }

    K0 = malloc(ncomp * sizeof(*K0));
    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);
    phase = flash_calculation_phase_new(eos, comp_X);

    for (i = 0; i < n_pres; i++) {
        K = NULL;
        Fv = Fv0 = 0.0;

        for (j = 0; j < n_temp; j++) {
            eos->pres = pres_list[i];
            eos->temp = temp_list[j];

            sm->pres[i * n_temp + j] = pres_list[i];
            sm->temp[i * n_temp + j] = temp_list[j];

            status = flash_calculation_stability_analysis_QNSS(phase, K0, 1e-10);

            if (status == 1) {
                if (phase->phase_no == 0) {
                    sm->F[i * n_temp + j] = 0.0;
                }
                else {
                    sm->F[i * n_temp + j] = 1.0;
                }

                for (k = 0; k < ncomp; k++) {
                    sm->K[i * n_temp + j][k] = 0.0;
                }
            }
            else if (status == 0) {
                double *K00;

                K00 = malloc(ncomp * sizeof(*K00));

                if (fsa == NULL) {
                    Fv0 = Fv;

                    for(k = 0; k < ncomp; k++) {
                        if (K == NULL) {
                            K00[k] = K0[k];
                        }
                        else {
                            K00[k] = K[k];
                        }
                    }
                }
                else {
                    int flag = 0;
                    double input[ncomp + 2], K000[ncomp], prediction_time;

                    for (k = 0; k < ncomp; k++) {
                        input[k] = comp_X[k];
                    }
                    input[ncomp] = temp_list[j];
                    input[ncomp + 1] = pres_list[i];

#if 0
                    printf("input:");
                    for (k = 0; k < ncomp + 2; k++) {
                        printf("%e ", input[k]);
                    }
                    printf("\n");
#endif

                    prediction_time = flash_calculation_get_time(NULL);
                    flag = flash_calculation_split_ann_predict(fsa, input, ncomp + 2, 
                            &Fv0, K00);
                    prediction_time = flash_calculation_get_time(NULL) - prediction_time;
                    split_pred_time += prediction_time;

                    if (verb) {
                        printf("Prediction time: %e\n", prediction_time);
                    }

                    if (!flag) {
                        Fv0 = -1.0;
                    }
                    else {
                        if (Fv0 > 1.0) {
                            Fv0 = 0.9;
                        }

                        if (Fv0 < 0.0) {
                            Fv0 = 0.1;
                        }

#if 0
                        printf("Fv0: %e\n", Fv0);
                        printf("K: ");
#endif
                        flash_calculation_estimate_K(eos, K000);

                        for (k = 0; k < ncomp; k++) {
                            //printf("%e ", K00[k]);

                            if (K00[k] < 0.0 || fabs(log(K00[k])) < 1e-4) {
                                K00[k] = K000[k];
                            }
                        }
                        //printf("\n");
                    }
                }

                solve_time = flash_calculation_get_time(NULL);
                Fv = flash_calculation_two_phase_flash_Calculation_QNSS(eos, 
                        comp_X, K00, Fv0, 1e-10);
                solve_time = flash_calculation_get_time(NULL) - solve_time;
                split_solve_time += solve_time;

                if (verb && ((fabs(Fv) < 1e-5 || fabs(Fv - 1.0) < 1e-5))) {
                    printf("Split calculation time: %e\n", solve_time);
                }

#if 0
                printf("FFv: %e\n", Fv);
                printf("FK: \n");
#endif
                for(k = 0; k < ncomp; k++) {
                    K0[k] = K00[k];

#if 0
                    printf("%e ", K00[k]);
#endif
                }
                free(K00);
                //printf("\n");

                sm->F[i * n_temp + j] = Fv;
                for (k = 0; k < ncomp; k++) {
                    sm->K[i * n_temp + j][k] = K0[k];
                }

                K = sm->K[i * n_temp + j];
            }

            for (k = 0; k < ncomp; k++) {
                sm->x[i * n_temp + j][k] = comp_X[k];
            }
        }
    }

    if (output_name != NULL) {
        flash_calculation_output_split_calculation_map(sm, comp_X, 
                ncomp, 0, output_name);
    }

    free(pres_list);
    free(temp_list);
    flash_calculation_phase_free(&phase);
    free(eos);

    return sm;
}

SPLIT_PM_MAP * flash_calculation_draw_split_calculation_map_PM(COMP_LIST *comp_list, 
        double *comp_X, double T, double P_min, double P_max, double dP, 
        int selected_component, double dxx, FLASH_SPLIT_ANN *fsa, char *output_name)
{
    int i, j, k, ncomp = comp_list->ncomp;
    int n_pres, n_x;
    double *pres_list, *x_list;
    SPLIT_PM_MAP *sm;
    int first, status;
    double Fv, Fv0, *K0, *K, *x, sum_no_selected, 
           solve_time;
    EOS *eos;
    PHASE *phase;
    char file_name[100];
    FILE *fp;

    n_pres = (int)((P_max - P_min) / dP);
    pres_list = malloc(n_pres * sizeof(*pres_list));
    for (i = 0; i < n_pres; i++) {
        pres_list[i] = P_min + i * dP;
    }

    n_x = (int)((1.0 - 0.001) / dxx) + 1;
    x_list = malloc(n_x * sizeof(*x_list));
    for (i = 0; i < n_x - 1; i++) {
        x_list[i] = 0.001 + i * dxx;
    }
    x_list[n_x - 1] = 0.999;

    sm = malloc(sizeof(*sm));
    sm->n = n_pres * n_x;
    sm->pres = malloc(n_pres * n_x* sizeof(*(sm->pres)));
    sm->F = malloc(n_pres * n_x* sizeof(*(sm->F)));
    sm->K = malloc(n_pres * n_x* sizeof(*(sm->K)));
    for (i = 0; i < n_pres * n_x; i++) {
        sm->K[i] = malloc(ncomp * sizeof(*(sm->K[i])));
    }
    sm->x = malloc(n_pres * n_x* sizeof(*(sm->x)));
    for (i = 0; i < n_pres * n_x; i++) {
        sm->x[i] = malloc(ncomp * sizeof(*(sm->x[i])));
    }

    K0 = malloc(ncomp * sizeof(*K0));
    x = malloc(ncomp * sizeof(*x));
    eos = flash_calculation_EOS_new(comp_list, 0.0, T, 0);
    phase = flash_calculation_phase_new(eos, x);

    sum_no_selected = 0.0;
    for (i = 0; i < ncomp; i++) {
        if (i != selected_component) {
            sum_no_selected += comp_X[i];
        }
    }

    for (i = 0; i < n_pres; i++) {
        K = NULL;
        Fv = Fv0 = 0.0;

        for (j = 0; j < n_x; j++) {
            for (k = 0; k < ncomp; k++) {
                if (k == selected_component) {
                    x[k] = x_list[j];
                }
                else {
                    x[k] = (1.0 - x_list[j]) * comp_X[k] / sum_no_selected;
                }
            }

            for (k = 0; k < ncomp; k++) {
                sm->x[i * n_x + j][k] = x[k];
            }

            eos->pres = pres_list[i];
            eos->temp = T;

            sm->pres[i * n_x + j] = pres_list[i];

            status = flash_calculation_stability_analysis_QNSS(phase, K0, 1e-10);

            if (status == 1) {
                if (phase->phase_no == 0) {
                    sm->F[i * n_x + j] = 0.0;
                }
                else {
                    sm->F[i * n_x + j] = 1.0;
                }

                for (k = 0; k < ncomp; k++) {
                    sm->K[i * n_x + j][k] = 0.0;
                }
            }
            else if (status == 0) {
                double *K00;

                K00 = malloc(ncomp * sizeof(*K00));

                if (fsa == NULL) {
                    Fv0 = Fv;

                    for(k = 0; k < ncomp; k++) {
                        if (K == NULL) {
                            K00[k] = K0[k];
                        }
                        else {
                            K00[k] = K[k];
                        }
                    }
                }
                else {
                    int flag = 0;
                    double input[ncomp + 1], K000[ncomp], prediction_time;

                    for (k = 0; k < ncomp; k++) {
                        input[k] = x[k];
                    }
                    input[ncomp] = pres_list[i];

                    prediction_time = flash_calculation_get_time(NULL);
                    flag = flash_calculation_split_ann_predict(fsa, input, ncomp + 1, 
                            &Fv0, K00);
                    prediction_time = flash_calculation_get_time(NULL) - prediction_time;
                    split_pred_time += prediction_time;

                    if (verb) {
                        printf("Prediction time: %e\n", prediction_time);
                    }

                    if (!flag) {
                        Fv0 = -1.0;
                    }
                    else {
                        if (Fv0 > 1.0) {
                            Fv0 = 0.9;
                        }

                        if (Fv0 < 0.0) {
                            Fv0 = 0.1;
                        }

                        if (verb) {
                        printf("Fv0: %e\n", Fv0);
                        printf("K: ");
                        }

                        flash_calculation_estimate_K(eos, K000);

                        for (k = 0; k < ncomp; k++) {
                            if (verb) {
                            printf("%e ", K00[k]);
                            }

                            if (K00[k] < 0.0 || fabs(log(K00[k])) < 1e-4) {
                                K00[k] = K000[k];
                            }
                        }

                        if (verb) {
                        printf("\n");
                        }
                    }
                }

                solve_time = flash_calculation_get_time(NULL);
                Fv = flash_calculation_two_phase_flash_Calculation_QNSS(eos, 
                        x, K00, Fv0, 1e-10);
                solve_time = flash_calculation_get_time(NULL) - solve_time;
                split_solve_time += solve_time;

                if (verb && ((fabs(Fv) < 1e-5 || fabs(Fv - 1.0) < 1e-5))) {
                    printf("Split calculation time: %e\n", solve_time);
                }


                if (verb) {
                printf("Fv: %e\n", Fv);
                printf("KK: ");
                }
                for(k = 0; k < ncomp; k++) {
                    K0[k] = K00[k];

                    if (verb) {
                    printf("%e ", K00[k]);
                    }
                }
                free(K00);
                if (verb) {
                printf("\n");
                }

                sm->F[i * n_x + j] = Fv;
                for (k = 0; k < ncomp; k++) {
                    sm->K[i * n_x + j][k] = K0[k];
                }

                K = sm->K[i * n_x + j];
            }
        }
    }

    if (output_name != NULL) {
        flash_calculation_output_split_calculation_map_PM(sm, NULL, 
                ncomp, 0, output_name);
    }

    free(pres_list);
    free(x_list);
    free(K0);
    flash_calculation_phase_free(&phase);
    free(eos);

    return sm;
}

void flash_calculation_phase_envelope_output(PHASE_ENVELOPE *pe,
        double *comp_X, int ncomp, char *output_name)
{
    char file_name[100], file_name_upper[100], file_name_down[100];
    FILE *fp, *fp_upper, *fp_down;
    double T, P;
    int i, j;
    int flag, flag1;

    if (output_name == NULL) 
        return;

    sprintf(file_name, "%s-phase-envelope.csv", output_name);
    fp = fopen(file_name, "a");

    for (j = 0; j < ncomp; j++) {
        fprintf(fp, "Component %d,", j);
    }
    fprintf(fp, "Temperature,Pressure\n");

    flag = 0; 
    flag1 = 0;
    for (i = 0; i < pe->n; i++) {
        if (pe->Ps[i] > 1.0) {
            flag = 1;
        }

        if (!flag && !flag1 && i % 5 != 0) {
            flag1 = flag;

            continue;
        }
        flag1 = flag;

        if (comp_X != NULL) {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", comp_X[j]);
            }
        }
        else {
            for (j = 0; j < ncomp; j++) {
                fprintf(fp, "%lf,", pe->xs[i][j]);
            }
        }

        fprintf(fp, "%lf,%lf\n", pe->Ts[i], pe->Ps[i]);

    }
    fclose(fp);

    sprintf(file_name_upper, "%s-phase-envelope-upper.csv", output_name);
    fp_upper = fopen(file_name_upper, "a");

    sprintf(file_name_down, "%s-phase-envelope-down.csv", output_name);
    fp_down = fopen(file_name_down, "a");

    flag = 0; 
    flag1 = 0;
    for (i = 0; i < pe->n; i++) {
        T = pe->Ts[i];
        P = pe->Ps[i];

        if (pe->Ps[i] > 1.0) {
            flag = 1;
        }
        else {
            flag = 0;
        }

        if (!flag && !flag1 && i % 5 != 0 
                && (i != pe->n - 1)) {
            flag1 = flag;

            continue;
        }
        flag1 = flag;

        if (i > 1) {
            if (fabs(T - pe->Ts[i - 1]) > 1e-5) {
                double dP, dT, rate, rate0;

                dP = P - pe->Ps[i - 1];
                dT = T - pe->Ts[i - 1];
                rate = fabs(dP / dT);

                dP = pe->Ps[i - 2] - pe->Ps[i - 1];
                dT = pe->Ts[i - 2] - pe->Ts[i - 1];

                if (dT == 0.0) {
                    rate0 = rate;
                }
                else {
                    rate0 = fabs(dP / dT);
                }

                if (rate > rate0 + 30.0 || rate + 30.0 < rate0) {
                    continue;
                }
            }
        }

        if (i == 0 || (i > 0 && T > pe->Ts[i - 1])) {
            if (comp_X != NULL) {
                for (j = 0; j < ncomp; j++) {
                    fprintf(fp_upper, "%lf,", comp_X[j]);
                }
            }
            else {
                for (j = 0; j < ncomp; j++) {
                    fprintf(fp_upper, "%lf,", pe->xs[i][j]);
                }
            }

            fprintf(fp_upper, "%lf,%lf\n", pe->Ts[i], pe->Ps[i]);
        }
        else {
            if (comp_X != NULL) {
                for (j = 0; j < ncomp; j++) {
                    fprintf(fp_down, "%lf,", comp_X[j]);
                }
            }
            else {
                for (j = 0; j < ncomp; j++) {
                    fprintf(fp_down, "%lf,", pe->xs[i][j]);
                }
            }

            fprintf(fp_down, "%lf,%lf\n", pe->Ts[i], pe->Ps[i]);
        }
    }

    fclose(fp_upper);
    fclose(fp_down);
}

/*
# ## 9. Draw Phase Diagram
# For a feed composition, this function will draw a phase diagram at the given 
vapour mole fraction, temperature range and pressure range. For the example, 
the following figure shows the phase diagram at $F_v = [0.3, 0.5, 0.8]$ for 
the temperature in [200 K, 340 K] and the pressure in [20 atm, 80 atm].
*/

PHASE_DIAGRAM * flash_calculation_draw_phase_diagram(COMP_LIST *comp_list, 
        double *comp_X, double T_min, double T_max, double P_min, double P_max, 
        double *Fv_list, int nF, double dT, double dP, PHASE_ENVELOPE *pe, 
        CRITICAL_POINT *cp, char *output_name)
{
    PHASE_DIAGRAM *pd;
    char file_name[100];
    FILE *fp;
    int i, j, k, ncomp = comp_list->ncomp;
    EOS *eos;

    eos = flash_calculation_EOS_new(comp_list, 0.0, 0.0, 0);

    pd = flash_calculation_phase_diagram_construction(eos, comp_X, 
            P_min, P_max, T_min, T_max, Fv_list, nF, pe, cp, 
            dP, dT);

    flash_calculation_phase_envelope_output(pd->pe,
            comp_X, ncomp, output_name);

    if (output_name == NULL) {
        for (j = 0; j < pd->n_line; j++) {
            sprintf(file_name, "%s-phase-diagram-F%lf.csv", output_name, pd->pl[j]->F);
            fp = fopen(file_name, "a");

            for (i = 0; i < ncomp; i++) {
                fprintf(fp, "Component %d,", i);
            }
            fprintf(fp, "Temperature,Pressure,Fv\n");

            for (i = 0; i < pd->pl[j]->n; i++) {
                for (k = 0; k < ncomp; k++) {
                    fprintf(fp, "%lf,", comp_X[k]);
                }

                fprintf(fp, "%lf,%lf,%lf\n", pd->pl[j]->T[i], pd->pl[j]->P[i], pd->pl[j]->F);
            }
            fclose(fp);
        }
    }

    printf("Critical Point: (%lf, %lf)\n", pd->cp->Tc, pd->cp->Pc);

    free(eos);

    return pd;
}

FLASH_MODEL * flash_calculation_model_new(char *model_file)
{
    FILE *fn;
    char line[10240];
    CSV_LINE **csv_value, *csv_line0; 
    int nline, i, j;
    FLASH_MODEL *fm;

    fm = malloc(sizeof(*fm));
    fm->type = malloc(100 * sizeof(*(fm->type)));
    fm->prop_file = malloc(100 * sizeof(*(fm->prop_file)));
    fm->bin_file = malloc(100 * sizeof(*(fm->bin_file)));
    fm->z_file = malloc(100 * sizeof(*(fm->z_file)));
    fm->output = malloc(100 * sizeof(*(fm->output)));
    fm->nF = 0;
    fm->F = NULL;
    fm->split_ann = NULL;
    fm->split_ann_level = 0;

    csv_value = malloc(20 * sizeof(*csv_value));

    fn = fopen(model_file, "r");
    nline = 0;
    while (fgets(line, 10240, fn)){
        csv_line0 = get_CSV_line(line, NULL);
        csv_value[nline] = csv_line0;

        nline++;
    }
    fclose(fn);

    for (i = 0; i < nline; i++) {
        char *name;

        name = csv_value[i]->value[0];

        if (strcmp(name, "type") == 0) {
            strcpy(fm->type, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "prop_file") == 0) {
            strcpy(fm->prop_file, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "bin_file") == 0) {
            strcpy(fm->bin_file, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "z_file") == 0) {
            strcpy(fm->z_file, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "T_min") == 0) {
            fm->T_min = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "T_max") == 0) {
            fm->T_max = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "P_min") == 0) {
            fm->P_min = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "P_max") == 0) {
            fm->P_max = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "dT") == 0) {
            fm->dT = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "dP") == 0) {
            fm->dP = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "dx") == 0) {
            fm->dx = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "F") == 0) {
            fm->nF = csv_value[i]->n - 1;
            fm->F = malloc(fm->nF * sizeof(*(fm->F)));

            for (j = 0; j < fm->nF; j++) {
                fm->F[j] = atof(csv_value[i]->value[j + 1]);
            }
        }
        else if (strcmp(name, "T") == 0) {
            fm->T = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "selected_component") == 0) {
            fm->selected_component = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "dxx") == 0) {
            fm->dxx = atof(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "split_ann") == 0) {
            fm->split_ann = malloc(100 * sizeof(*(fm->split_ann)));

            strcpy(fm->split_ann, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "split_ann_level") == 0) {
            fm->split_ann_level = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "stab_ann") == 0) {
            fm->stab_ann = malloc(100 * sizeof(*(fm->stab_ann)));

            strcpy(fm->stab_ann, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "stab_ann_level") == 0) {
            fm->stab_ann_level = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "stab_ann_safeguard") == 0) {
            fm->stab_ann_safeguard = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "stab_ann_delta_p") == 0) {
            fm->stab_ann_delta_p = atoi(csv_value[i]->value[1]);
        }
        else if (strcmp(name, "output") == 0) {
            strcpy(fm->output, csv_value[i]->value[1]);
        }
        else if (strcmp(name, "verb") == 0) {
            fm->verb = atoi(csv_value[i]->value[1]);
        }
    }

    return fm;
}

void flash_calculation_phase_envelope_free(PHASE_ENVELOPE **pe)
{
    PHASE_ENVELOPE *pe0;

    pe0 = *pe;
    free(pe0->Ps);
    free(pe0->Ts);

    free(*pe);
}

void flash_calculation_phase_envelope_pm_free(PHASE_ENVELOPE_PM **pe_pm)
{
    int i;
    PHASE_ENVELOPE_PM *pe_pm0 = *pe_pm;

    free(pe_pm0->Ps);

    for (i = 0; i < pe_pm0->n; i++) {
        free(pe_pm0->xs[i]);
    }
    free(pe_pm0->xs);

    free(*pe_pm);
}

void flash_calculation_critical_point_free(CRITICAL_POINT **cp)
{
    free(*cp);
}

void flash_calculation_phase_diagram_free(PHASE_DIAGRAM **pd)
{
    int i;
    PHASE_DIAGRAM *pd0 = *pd;

    flash_calculation_phase_envelope_free(&(pd0->pe));
    flash_calculation_critical_point_free(&(pd0->cp));

    for (i = 0; i < pd0->n_line; i++) {
        free(pd0->pl[i]->P);
        free(pd0->pl[i]->T);
        free(pd0->pl[i]);
    }

    free(pd0->pl);

    free(*pd);
}

void flash_calculation_stability_map_free(STABILITY_MAP **sm)
{
    int i;
    STABILITY_MAP *sm0 = *sm;

    free(sm0->unstable_pres);
    free(sm0->unstable_temp);

    free(sm0->liquid_pres);
    free(sm0->liquid_temp);

    free(sm0->vapor_pres);
    free(sm0->vapor_temp);

    free(*sm);
}

void flash_calculation_stability_PM_map_free(STABILITY_PM_MAP **sm)
{
    int i;
    STABILITY_PM_MAP *sm0 = *sm;

    for (i = 0; i < sm0->n_unstable; i++) {
        free(sm0->unstable_x[i]);
    }
    free(sm0->unstable_pres);
    free(sm0->unstable_x);

    for (i = 0; i < sm0->n_unstable; i++) {
        free(sm0->liquid_x[i]);
    }
    free(sm0->liquid_pres);
    free(sm0->liquid_x);

    for (i = 0; i < sm0->n_unstable; i++) {
        free(sm0->vapor_x[i]);
    }
    free(sm0->vapor_pres);
    free(sm0->vapor_x);

    free(*sm);
}

void flash_calculation_split_map_free(SPLIT_MAP **sm)
{
    int i;
    SPLIT_MAP *sm0 = *sm;

    free(sm0->temp);
    free(sm0->pres);
    free(sm0->F);

    for (i = 0; i < sm0->n; i++) {
        free(sm0->K[i]);
        free(sm0->x[i]);
    }

    free(sm0->K);

    free(*sm);
}

void flash_calculation_split_PM_map_free(SPLIT_PM_MAP **sm)
{
    int i;
    SPLIT_PM_MAP *sm0 = *sm;

    free(sm0->pres);
    free(sm0->F);

    for (i = 0; i < sm0->n; i++) {
        free(sm0->K[i]);
        free(sm0->x[i]);
    }

    free(sm0->K);
    free(sm0->x);

    free(*sm);
}

void flash_calculation_flash_model_free(FLASH_MODEL **fm)
{
    FLASH_MODEL *fm0 = *fm;

    free(fm0->type);
    free(fm0->prop_file);
    free(fm0->bin_file);
    free(fm0->z_file);
    if (fm0->F != NULL) {
        free(fm0->F);
    }
    free(fm0->output);

    if (fm0->split_ann != NULL)
        free(fm0->split_ann);

    free(*fm);
}

static int flash_calculation_generate_x_number(double start,
        double end, double dx, int ncomp)
{
    int i, sum, sum2;

    if (ncomp == 2) {
        sum = (int)((end - start) / dx);

        return sum + 1;
    }

    sum2 = (int)((end - start) / dx);
    sum = 0;
    for (i = 0; i < sum2; i++) {
        int sum3;
        double x0;

        x0 = start + i * dx;
        sum3 = flash_calculation_generate_x_number(start, 
                end - x0, dx, ncomp - 1);

        sum += sum3;
    }

    sum += ncomp - 1;

    return sum;
}

static int flash_calculation_generate_x(double start, double end, double dx, 
        int ncomp, double min_x, double ***x_list)
{
    int i, j, k, count, sum, sum2, ncomp2;
    double **x_list0;

    if (ncomp == 2) {
        sum = (int)((end - start) / dx);

        *x_list = malloc((sum + 1)* sizeof(**x_list));
        x_list0 = *x_list;

        //printf("N222222: %lf, %lf\n", start, end);
        for (i = 0; i < sum; i++) {
            x_list0[i] = malloc(2 * sizeof(double));

            x_list0[i][0] = start + i * dx;
            //printf("[0]: %lf\n", x_list0[i][0]);
            x_list0[i][1] = end - (start + i * dx);
            //printf("[1]: %lf\n", x_list0[i][1]);
        }

        x_list0[sum] = malloc(2 * sizeof(double));
        x_list0[sum][1] = dx;
        x_list0[sum][0] = end - dx;

        return sum + 1;
    }

    sum = flash_calculation_generate_x_number(start,
            end, dx, ncomp);
    sum2 = (int)((end - start) / dx);

    //printf("number: %d\n", sum);

    *x_list = malloc(sum * sizeof(**x_list));
    x_list0 = *x_list;
    for (i = 0; i < sum; i++) {
        x_list0[i] = malloc(ncomp * sizeof(double));
    }

    //printf("ncomp: %d\n", ncomp);
    //printf("end: %e, start: %e, dx: %e\n", end, start, dx);
    //printf("sum2: %d\n", sum2);

    /* add list */
    count = 0;
    for (i = 0; i < sum2; i++) {
        int sum3;
        double x0, **x_list_b;

        x0 = start + i * dx;
        //printf("strat: %lf, x0: %lf, end: %lf, end-x0: %lf\n", 
        //        start, x0, end, end-x0);
        sum3 = flash_calculation_generate_x(start, end - x0, dx, 
                ncomp - 1, min_x, &x_list_b);
        //printf("sum3: %d\n", sum3);

        for (j = 0; j < sum3; j++) {
            x_list0[count][0] = x0;

            for (k = 0; k < ncomp - 1; k++) {
                x_list0[count][k + 1] = x_list_b[j][k];
            }

            count += 1;
            free(x_list_b[j]);
        }

        free(x_list_b);
    }

    while(1) {
        if (end - min_x * (ncomp - 1) > min_x) {
            break;
        }
        else {
            min_x *= 0.5;
        }
    }

    for (i = 0; i < ncomp - 1; i++) {
        for (j = 0; j < ncomp; j++) {
            if (i == j) {
                x_list0[count][j] = end - min_x * (ncomp - 1);
            }
            else {
                x_list0[count][j] = min_x;
            }
        }
        count += 1;
    }

    if (count != sum)
        printf("JAKDFJSLKDJFSLKDJFSDLKFJLSDKJFSDKJFSLDKJFSDLKFJD\n");

    return sum;
}

void flash_calculation_generate_stability_analysis_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, FLASH_STAB_ANN *fsa, char *output)
{
    int i, j, ncomp = comp_list->ncomp;
    STABILITY_MAP *sm;
    double *x;

    x = malloc(ncomp * sizeof(*x));

    printf("   -- Generating stability analysis data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);
        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        sm = flash_calculation_draw_stability_analysis_map(comp_list, x, 
                T_min, T_max, P_min, P_max, dT, dP, fsa, output);

        flash_calculation_stability_map_free(&sm);

        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");
    }
    printf("   -- Done\n");

    free(x);
}

void flash_calculation_generate_stability_analysis_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double P_min, double P_max, double dP, 
        double dxx, FLASH_STAB_ANN *fsa, char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    STABILITY_PM_MAP *sm0, *sm;
    double *x;
    double max_P = 0.0, min_P = 1e10;
    double max_P0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));

    sm = malloc(sizeof(*sm));
    sm->n_unstable = 0;
    sm->unstable_pres = malloc(sizeof(*(sm->unstable_pres)));
    sm->unstable_x = malloc(sizeof(*(sm->unstable_x)));
    sm->n_liquid = 0;
    sm->liquid_pres = malloc(sizeof(*(sm->liquid_pres)));
    sm->liquid_x = malloc(sizeof(*(sm->liquid_x)));
    sm->n_vapor = 0;
    sm->vapor_pres = malloc(sizeof(*(sm->vapor_pres)));
    sm->vapor_x = malloc(sizeof(*(sm->vapor_x)));

    printf("   -- Generating stability analysis data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);
        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        sm0 = flash_calculation_draw_stability_analysis_map_PM(comp_list, x, 
                T, P_min, P_max, dP, 0, dxx, fsa, output);

        sm->unstable_pres = realloc(sm->unstable_pres,
                (sm->n_unstable + sm0->n_unstable) * sizeof(*(sm->unstable_pres)));
        sm->liquid_pres = realloc(sm->liquid_pres,
                (sm->n_liquid + sm0->n_liquid) * sizeof(*(sm->liquid_pres)));
        sm->vapor_pres = realloc(sm->vapor_pres,
                (sm->n_vapor + sm0->n_vapor) * sizeof(*(sm->vapor_pres)));
        sm->unstable_x = realloc(sm->unstable_x,
                (sm->n_unstable + sm0->n_unstable) * sizeof(*(sm->unstable_x)));
        sm->liquid_x = realloc(sm->liquid_x,
                (sm->n_liquid + sm0->n_liquid) * sizeof(*(sm->liquid_x)));
        sm->vapor_x = realloc(sm->vapor_x,
                (sm->n_vapor + sm0->n_vapor) * sizeof(*(sm->vapor_x)));

        for (j = sm->n_unstable; j < sm->n_unstable + sm0->n_unstable; j++) {
            sm->unstable_pres[j] = sm0->unstable_pres[j - sm->n_unstable];

            sm->unstable_x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->unstable_x[j][k] = sm0->unstable_x[j - sm->n_unstable][k];
            }

            if (max_P < sm->unstable_pres[j]) {
                max_P = sm->unstable_pres[j];
            }

            if (min_P > sm->unstable_pres[j]) {
                min_P = sm->unstable_pres[j];
            }
        }
        sm->n_unstable += sm0->n_unstable;

        for (j = sm->n_liquid; j < sm->n_liquid + sm0->n_liquid; j++) {
            sm->liquid_pres[j] = sm0->liquid_pres[j - sm->n_liquid];

            sm->liquid_x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->liquid_x[j][k] = sm0->liquid_x[j - sm->n_liquid][k];
            }

            if (max_P < sm->liquid_pres[j]) {
                max_P = sm->liquid_pres[j];
            }

            if (min_P > sm->liquid_pres[j]) {
                min_P = sm->liquid_pres[j];
            }
        }
        sm->n_liquid += sm0->n_liquid;

        for (j = sm->n_vapor; j < sm->n_vapor + sm0->n_vapor; j++) {
            sm->vapor_pres[j] = sm0->vapor_pres[j - sm->n_vapor];

            sm->vapor_x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->vapor_x[j][k] = sm0->vapor_x[j - sm->n_vapor][k];
            }

            if (max_P < sm->vapor_pres[j]) {
                max_P = sm->vapor_pres[j];
            }

            if (min_P > sm->vapor_pres[j]) {
                min_P = sm->vapor_pres[j];
            }
        }
        sm->n_vapor += sm0->n_vapor;

        flash_calculation_stability_PM_map_free(&sm0);

        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");
    }
    printf("   -- Done\n");

    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
    }

    for (i = 0; i < sm->n_unstable; i++) {
        sm->unstable_pres[i] = (sm->unstable_pres[i] - min_P0) / (max_P0 - min_P0);
    }
    for (i = 0; i < sm->n_liquid; i++) {
        sm->liquid_pres[i] = (sm->liquid_pres[i] - min_P0) / (max_P0 - min_P0);
    }
    for (i = 0; i < sm->n_vapor; i++) {
        sm->vapor_pres[i] = (sm->vapor_pres[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_output_stability_analysis_map_PM(sm, 
            NULL, ncomp, file_scaled);

    free(sm->unstable_pres);
    for (i = 0; i < sm->n_unstable; i++) {
        free(sm->unstable_x[i]);
    }
    free(sm->unstable_x);

    free(sm->liquid_pres);
    for (i = 0; i < sm->n_liquid; i++) {
        free(sm->liquid_x[i]);
    }
    free(sm->liquid_x);

    free(sm->vapor_pres);
    for (i = 0; i < sm->n_vapor; i++) {
        free(sm->vapor_x[i]);
    }
    free(sm->vapor_x);

    free(sm);
    free(x);
}

void flash_calculation_generate_split_calculation_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, FLASH_SPLIT_ANN *fsa, char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    SPLIT_MAP *sm0, *sm;
    double *x;
    double max_T = 0.0, max_P = 0.0, 
           min_T = 1e10, min_P = 1e10;
    double max_T0, max_P0, min_T0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));

    sm = malloc(sizeof(*sm));
    sm->n = 0;
    sm->temp = malloc(sizeof(*(sm->temp)));
    sm->pres = malloc(sizeof(*(sm->pres)));
    sm->F = malloc(sizeof(*(sm->F)));
    sm->K = malloc(sizeof(*(sm->K)));
    sm->x = malloc(sizeof(*(sm->x)));

    printf("   -- Generating split calculation data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);

        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        sm0 = flash_calculation_draw_split_calculation_map(comp_list, x, 
                T_min, T_max, P_min, P_max, dT, dP, fsa, output);

        sm->temp = realloc(sm->temp, 
                (sm->n + sm0->n) * sizeof(*(sm->temp)));
        sm->pres = realloc(sm->pres, 
                (sm->n + sm0->n) * sizeof(*(sm->pres)));
        sm->F = realloc(sm->F, 
                (sm->n + sm0->n) * sizeof(*(sm->F)));
        sm->K = realloc(sm->K, 
                (sm->n + sm0->n) * sizeof(*(sm->K)));
        sm->x = realloc(sm->x, 
                (sm->n + sm0->n) * sizeof(*(sm->x)));

        for (j = sm->n; j < sm->n + sm0->n; j++) {
            sm->temp[j] = sm0->temp[j - sm->n];
            sm->pres[j] = sm0->pres[j - sm->n];
            sm->F[j] = sm0->F[j - sm->n];

            sm->K[j] = malloc(ncomp * sizeof(double));
            sm->x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->K[j][k] = sm0->K[j - sm->n][k];
                sm->x[j][k] = x[k];
            }

            if (max_T < sm->temp[j]) {
                max_T = sm->temp[j];
            }
            if (max_P < sm->pres[j]) {
                max_P = sm->pres[j];
            }

            if (min_T > sm->temp[j]) {
                min_T = sm->temp[j];
            }
            if (min_P > sm->pres[j]) {
                min_P = sm->pres[j];
            }
        }
        sm->n += sm0->n;

        flash_calculation_split_map_free(&sm0);
        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");

    }
    printf("   -- Done\n");

    sprintf(file_scaled, "%s-non-scaled", output);
    flash_calculation_output_split_calculation_map(sm, NULL, 
            ncomp, 1, file_scaled);

    MPI_Allreduce(&max_T, &max_T0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_T, &min_T0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
        printf("Minimal Temperature: %e, Maximal Temperature: %e\n", 
                min_T0, max_T0);
    }

    for (i = 0; i < sm->n; i++) {
        sm->temp[i] = (sm->temp[i] - min_T0) / (max_T0 - min_T0);
        sm->pres[i] = (sm->pres[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_output_split_calculation_map(sm, NULL, 
            ncomp, 1, file_scaled);

    free(sm->pres);
    free(sm->temp);
    free(sm->F);
    for (i = 0; i < sm->n; i++) {
        free(sm->x[i]);
        free(sm->K[i]);
    }
    free(sm->x);
    free(sm->K);
    free(sm);

    free(x);
}

void flash_calculation_generate_split_calculation_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double P_min, double P_max, double dP, 
        double dxx, FLASH_SPLIT_ANN *fsa, char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    SPLIT_PM_MAP *sm0, *sm;
    double *x;
    double max_P = 0.0, min_P = 1e10;
    double max_P0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));

    sm = malloc(sizeof(*sm));
    sm->n = 0;
    sm->pres = malloc(sizeof(*(sm->pres)));
    sm->F = malloc(sizeof(*(sm->F)));
    sm->K = malloc(sizeof(*(sm->K)));
    sm->x = malloc(sizeof(*(sm->x)));

    printf("   -- Generating split calculation data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);

        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        sm0 = flash_calculation_draw_split_calculation_map_PM(comp_list, x, 
                T, P_min, P_max, dP, 0, dxx, fsa, output);

        sm->pres = realloc(sm->pres, 
                (sm->n + sm0->n) * sizeof(*(sm->pres)));
        sm->F = realloc(sm->F, 
                (sm->n + sm0->n) * sizeof(*(sm->F)));
        sm->K = realloc(sm->K, 
                (sm->n + sm0->n) * sizeof(*(sm->K)));
        sm->x = realloc(sm->x, 
                (sm->n + sm0->n) * sizeof(*(sm->x)));

        for (j = sm->n; j < sm->n + sm0->n; j++) {
            sm->pres[j] = sm0->pres[j - sm->n];
            sm->F[j] = sm0->F[j - sm->n];

            sm->K[j] = malloc(ncomp * sizeof(double));
            sm->x[j] = malloc(ncomp * sizeof(double));
            for (k = 0; k < ncomp; k++) {
                sm->K[j][k] = sm0->K[j - sm->n][k];
                sm->x[j][k] = sm0->x[j - sm->n][k];
            }

            if (max_P < sm->pres[j]) {
                max_P = sm->pres[j];
            }

            if (min_P > sm->pres[j]) {
                min_P = sm->pres[j];
            }
        }
        sm->n += sm0->n;

        flash_calculation_split_PM_map_free(&sm0);
        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");

    }
    printf("   -- Done\n");

    sprintf(file_scaled, "%s-non-scaled", output);
    flash_calculation_output_split_calculation_map_PM(sm, NULL, 
            ncomp, 1, file_scaled);

    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
    }

    for (i = 0; i < sm->n; i++) {
        sm->pres[i] = (sm->pres[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_output_split_calculation_map_PM(sm, NULL, 
            ncomp, 1, file_scaled);

    free(sm->pres);
    free(sm->F);
    for (i = 0; i < sm->n; i++) {
        free(sm->x[i]);
        free(sm->K[i]);
    }
    free(sm->x);
    free(sm->K);
    free(sm);

    free(x);
}

void flash_calculation_generate_phase_envelope_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T_min, double T_max, double P_min, 
        double P_max, double dT, double dP, char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    PHASE_DIAGRAM *pd;
    double *x;
    CRITICAL_POINT *cp;
    PHASE_ENVELOPE *pe;
    double max_T = 0.0, max_P = 0.0, 
           min_T = 1e10, min_P = 1e10;
    double max_T0, max_P0, min_T0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));
    pe = malloc(sizeof(*pe));
    pe->n = 0;
    pe->Ps = malloc(sizeof(*(pe->Ps)));
    pe->Ts = malloc(sizeof(*(pe->Ts)));
    pe->xs = malloc(sizeof(*(pe->xs)));

    printf("   -- Generating phase envelope data:\n");

    for (i = 0; i < nx; i++) {
        cp = malloc(sizeof(*cp));
        cp->Tc = 0.0;
        cp->Pc = 0.0;

        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);
        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        pd = flash_calculation_draw_phase_diagram(comp_list, x,
                T_min, T_max, P_min, P_max, NULL, 0, dT, dP, 
                NULL, cp, output);

        pe->Ps = realloc(pe->Ps, (pe->n + pd->pe->n) * sizeof(*(pe->Ps)));
        pe->Ts = realloc(pe->Ts, (pe->n + pd->pe->n) * sizeof(*(pe->Ts)));
        pe->xs = realloc(pe->xs, (pe->n + pd->pe->n) * sizeof(*(pe->xs)));

        for (j = pe->n; j < pe->n + pd->pe->n; j++) {
            pe->Ps[j] = pd->pe->Ps[j - pe->n];
            pe->Ts[j] = pd->pe->Ts[j - pe->n];
            pe->xs[j] = malloc(ncomp * sizeof(*(pe->xs[j])));
            for (k = 0; k < ncomp; k++) {
                pe->xs[j][k] = x[k];
            }

            if (max_T < pe->Ts[j]) {
                max_T = pe->Ts[j];
            }
            if (max_P < pe->Ps[j]) {
                max_P = pe->Ps[j];
            }

            if (min_T > pe->Ts[j]) {
                min_T = pe->Ts[j];
            }
            if (min_P > pe->Ps[j]) {
                min_P = pe->Ps[j];
            }
        }
        pe->n += pd->pe->n;

        flash_calculation_phase_diagram_free(&pd);

        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");
    }
    printf("   -- Done\n");

    MPI_Allreduce(&max_T, &max_T0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_T, &min_T0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
        printf("Minimal Temperature: %e, Maximal Temperature: %e\n", 
                min_T0, max_T0);
    }

    for (i = 0; i < pe->n; i++) {
        pe->Ts[i] = (pe->Ts[i] - min_T0) / (max_T0 - min_T0);
        pe->Ps[i] = (pe->Ps[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_phase_envelope_output(pe, NULL, ncomp, file_scaled);

    free(pe->Ts);
    free(pe->Ps);
    for (i = 0; i < pe->n; i++) {
        free(pe->xs[i]);
    }
    free(pe->xs);
    free(pe);

    free(x);
}

void flash_calculation_generate_phase_envelope_PM_data(COMP_LIST *comp_list, 
        int nx, double **x_list, double T, double dP, double dxx, 
        char *output)
{
    int i, j, k, ncomp = comp_list->ncomp;
    PHASE_ENVELOPE_PM *pe_pm0, *pe_pm;
    double *x;
    double max_P = 0.0, min_P = 1e10, max_P0, min_P0;
    char file_scaled[100];
    int myrank;

    x = malloc(ncomp * sizeof(*x));

    pe_pm = malloc(sizeof(*pe_pm));
    pe_pm->n = 0;
    pe_pm->Ps = malloc(sizeof(*(pe_pm->Ps)));
    pe_pm->xs = malloc(sizeof(*(pe_pm->xs)));

    printf("   -- Generating phase envelope data:\n");

    for (i = 0; i < nx; i++) {
        for (j = 0; j < ncomp; j++) {
            x[j] = x_list[i][j];
        }

        printf("      -- Data Set %d / %d\n", i, nx);
        printf("      --    x:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("\n");

        //for (j = 0; j < ncomp; j++) {
        pe_pm0 = flash_calculation_phase_saturation_envelope_construction_PM(comp_list, 
                x, T, 100.0, dP, 0, dxx, output);

        pe_pm->Ps = realloc(pe_pm->Ps, 
                (pe_pm->n + pe_pm0->n) * sizeof(*(pe_pm->Ps)));
        pe_pm->xs = realloc(pe_pm->xs, 
                (pe_pm->n + pe_pm0->n) * sizeof(*(pe_pm->xs)));

        for (j = pe_pm->n; j < pe_pm->n + pe_pm0->n; j++) {
            pe_pm->Ps[j] = pe_pm0->Ps[j - pe_pm->n];
            pe_pm->xs[j] = malloc(ncomp * sizeof(*(pe_pm->xs[j])));
            for (k = 0; k < ncomp; k++) {
                pe_pm->xs[j][k] = pe_pm0->xs[j - pe_pm->n][k];
            }

            if (max_P < pe_pm->Ps[j]) {
                max_P = pe_pm->Ps[j];
            }
            if (min_P > pe_pm->Ps[j]) {
                min_P = pe_pm->Ps[j];
            }
        }
        pe_pm->n += pe_pm0->n;

        flash_calculation_phase_envelope_pm_free(&pe_pm0);
        //}

        printf("      --    X:");
        for (j = 0; j < ncomp; j++) {
            printf("%lf ", x[j]);
        }
        printf("      --    Done\n");
    }
    printf("   -- Done\n");

    MPI_Allreduce(&max_P, &max_P0, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&min_P, &min_P0, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    if (myrank == 0) {
        printf("Minimal Pressure: %e, Maximal Pressure: %e\n", 
                min_P0, max_P0);
    }

    for (i = 0; i < pe_pm->n; i++) {
        pe_pm->Ps[i] = (pe_pm->Ps[i] - min_P0) / (max_P0 - min_P0);
    }

    sprintf(file_scaled, "%s-scaled", output);
    flash_calculation_phase_envelope_PM_output(pe_pm, ncomp, 0, file_scaled);

    free(pe_pm->Ps);
    for (i = 0; i < pe_pm->n; i++) {
        free(pe_pm->xs[i]);
    }
    free(pe_pm->xs);
    free(pe_pm);

    free(x);
}

int main(int argc, char **argv) {
    char *model_file, *prop_file, *binary_file, *z_file;
    COMP_LIST *comp_list;
    double pres, temp, *x, *K;
    int i, j;
    int nprocs, myrank;
    double Fv_list[3] = {0.2, 0.5, 0.8};
    FLASH_MODEL *fm;
    FLASH_SPLIT_ANN *fsa = NULL;
    FLASH_STAB_ANN *fsta = NULL;
    double mem_peak;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    printf("RANK[%d]: argc: %d, argv[1]: %s\n", myrank, argc, argv[1]);

    model_file = argv[1];
    fm = flash_calculation_model_new(model_file);

    verb = fm->verb;

    prop_file = fm->prop_file;
    binary_file = fm->bin_file;
    z_file = fm->z_file;

    printf("prop_name: %s\n", prop_file);
    printf("binary_name: %s\n", binary_file);
    printf("z_name: %s\n", z_file);

    comp_list = flash_calculation_component_new(prop_file, binary_file);
    x = flash_calculation_composition_new(z_file);

    if (fm->split_ann != NULL && fm->split_ann_level >= 0) {
        fsa = flash_calculation_split_ann_model_new(fm->split_ann, 
                fm->split_ann_level, comp_list->ncomp);
    }

    if (fm->stab_ann != NULL && fm->stab_ann_level >= 0) {
        fsta = flash_calculation_stab_ann_model_new(fm->stab_ann, 
                fm->stab_ann_level, fm->stab_ann_safeguard, 
                fm->stab_ann_delta_p);
    }

    if (strcmp(fm->type, "stability") == 0) {
        STABILITY_MAP *sm;

        sm = flash_calculation_draw_stability_analysis_map(comp_list, x, 
                fm->T_min, fm->T_max, fm->P_min, fm->P_max, fm->dT, fm->dP, 
                fsta, fm->output);

        flash_calculation_stability_map_free(&sm);
    }

    if (strcmp(fm->type, "stability_PM") == 0) {
        STABILITY_PM_MAP *sm;

        sm = flash_calculation_draw_stability_analysis_map_PM(comp_list, x, 
                fm->T, fm->P_min, fm->P_max, fm->dP, fm->selected_component, 
                fm->dxx, fsta, fm->output);

        flash_calculation_stability_PM_map_free(&sm);
    }

    if (strcmp(fm->type, "stability_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);

        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%04d", fm->output, myrank);
        }

        flash_calculation_generate_stability_analysis_data(comp_list, nx_rank, 
                x_list + nx_begin, fm->T_min, fm->T_max, fm->P_min, fm->P_max, 
                fm->dT, fm->dP, fsta, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }

    if (strcmp(fm->type, "stability_PM_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);

        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%04d", fm->output, myrank);
        }

        flash_calculation_generate_stability_analysis_PM_data(comp_list, nx_rank, 
                x_list + nx_begin, fm->T, fm->P_min, fm->P_max, fm->dP, 
                fm->dxx, fsta, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }

    if (strcmp(fm->type, "split") == 0) {
        SPLIT_MAP *sm;
        char file_scaled[100];

        sm = flash_calculation_draw_split_calculation_map(comp_list, x, 
                fm->T_min, fm->T_max, fm->P_min, fm->P_max, fm->dT, fm->dP, 
                fsa, fm->output);

        for (i = 0; i < sm->n; i++) {
            sm->temp[i] = (sm->temp[i] - fm->T_min) / (fm->T_max - fm->dT - fm->T_min);
            sm->pres[i] = (sm->pres[i] - fm->P_min) / (fm->P_max - fm->dP - fm->P_min);
        }

        sprintf(file_scaled, "%s-scaled", fm->output);
        flash_calculation_output_split_calculation_map(sm, x, 
                comp_list->ncomp, 1, file_scaled);

        flash_calculation_split_map_free(&sm);
    }

    if (strcmp(fm->type, "split_PM") == 0) {
        SPLIT_PM_MAP *sm;
        char file_scaled[100];

        sm = flash_calculation_draw_split_calculation_map_PM(comp_list, x, 
                fm->T, fm->P_min, fm->P_max, fm->dP, fm->selected_component,
                fm->dxx, fsa, fm->output);

        for (i = 0; i < sm->n; i++) {
            sm->pres[i] = (sm->pres[i] - fm->P_min) / (fm->P_max - fm->dP - fm->P_min);
        }

        sprintf(file_scaled, "%s-scaled", fm->output);
        flash_calculation_output_split_calculation_map_PM(sm, NULL, 
                comp_list->ncomp, 1, file_scaled);

        flash_calculation_split_PM_map_free(&sm);
    }

    if (strcmp(fm->type, "split_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);

        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%04d", fm->output, myrank);
        }

        flash_calculation_generate_split_calculation_data(comp_list, 
                nx_rank, x_list + nx_begin, fm->T_min, fm->T_max, fm->P_min, 
                fm->P_max, fm->dT, fm->dP, fsa, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }

    if (strcmp(fm->type, "split_PM_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);

        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%04d", fm->output, myrank);
        }

        flash_calculation_generate_split_calculation_PM_data(comp_list, 
                nx_rank, x_list + nx_begin, fm->T, fm->P_min, fm->P_max, 
                fm->dP, fm->dxx, fsa, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }


    if (strcmp(fm->type, "diagram") == 0) {
        PHASE_DIAGRAM *pd;

        pd = flash_calculation_draw_phase_diagram(comp_list, x, fm->T_min, fm->T_max,
                fm->P_min, fm->P_max, fm->F, fm->nF, fm->dT, fm->dP, NULL, NULL,
                fm->output);

        flash_calculation_phase_diagram_free(&pd);
    }

    if (strcmp(fm->type, "envelope_PM") == 0) {
        PHASE_ENVELOPE_PM *pe_pm;

        pe_pm = flash_calculation_phase_saturation_envelope_construction_PM(comp_list, 
                x, fm->T, 100.0, fm->dP, fm->selected_component, fm->dxx, fm->output);

        flash_calculation_phase_envelope_pm_free(&pe_pm);
    }

    if (strcmp(fm->type, "envelope_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);

        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        printf("rank[%d]---  nx: %d, nx_rank: %d, nx_begin: %d\n", myrank, nx, nx_rank, nx_begin);

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%03d", fm->output, myrank);
        }

        flash_calculation_generate_phase_envelope_data(comp_list, 
                nx_rank, x_list + nx_begin, fm->T_min, fm->T_max, fm->P_min, 
                fm->P_max, fm->dT, fm->dP, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }

    if (strcmp(fm->type, "envelope_PM_data") == 0) {
        int nx, nx_rank, nx_begin;
        double **x_list;
        char output_rank[100];

        nx = flash_calculation_generate_x(0.005, 1.0, fm->dx, comp_list->ncomp, 
                0.005, &x_list);

        nx_rank = nx / nprocs;
        nx_begin = nx_rank * myrank;
        if (myrank == nprocs - 1) {
            nx_rank += nx - nx_rank * nprocs;
        }

        printf("rank[%d]---  nx: %d, nx_rank: %d, nx_begin: %d\n", myrank, nx, nx_rank, nx_begin);

        if (nprocs == 1) {
            sprintf(output_rank, "%s", fm->output);
        }
        else {
            sprintf(output_rank, "%s-rank-%03d", fm->output, myrank);
        }

        flash_calculation_generate_phase_envelope_PM_data(comp_list, 
                nx_rank, x_list + nx_begin, fm->T, fm->dP, 
                fm->dxx, output_rank);

        for (i = 0; i < nx; i++) {
            free(x_list[i]);
        }
        free(x_list);
    }

    printf("Stability calculation iteration: %d\n", stab_itr);
    printf("Stability calculation ANN prediction time: %lf\n", stab_pred_time);
    printf("Stability calculation solve time: %lf\n", stab_solve_time);

    printf("Split calculation failure: %d\n", split_failure);
    printf("Split calculation iteration: %d\n", split_itr);
    printf("Split calculation ANN prediction time: %lf\n", split_pred_time);
    printf("Split calculation solve time: %lf\n", split_solve_time);

    if (fsa != NULL) {
        flash_calculation_split_ann_model_free(&fsa);
    }

    free(x);
    flash_calculation_flash_model_free(&fm);
    flash_calculation_component_free(&comp_list);

    //MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();

    return 0;
}



