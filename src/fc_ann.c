#include "fc.h"

static double split_pred_time = 0.;
static double stab_pred_time = 0.;

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
#if 0
        printf("x * W[%d]:\n", i);
        flash_calculation_tensor_dump(y[i]);
#endif

        flash_calculation_tensor_add(y[i], b[i], y[i]);
#if 0
        printf("x * W[%d] + b[%d]:\n", i, i);
        flash_calculation_tensor_dump(y[i]);
#endif

        flash_calculation_tensor_softmax(y[i], y[i]);
#if 0
        printf("softmax(x * W[%d] + b[%d]):\n", i, i);
        flash_calculation_tensor_dump(y[i]);
#endif

        x = y[i];
    }

    y[level - 1] = flash_calculation_tensor_matmul(x, W[level - 1], NULL);
#if 0
    printf("x * W[%d]:\n", level - 1, level - 1);
    flash_calculation_tensor_dump(y[level - 1]);
#endif

    flash_calculation_tensor_add(y[level - 1], b[level - 1], y[level - 1]);
#if 0
    printf("x * W[%d] + b[%d]:\n", level - 1, level - 1);
    flash_calculation_tensor_dump(y[level - 1]);
#endif

    value = y[level - 1]->value[0];
    //printf("value: %e\n", value);

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
    double pred_time;

    pred_time = flash_calculation_get_time(NULL);

    if (fsa->ann_F == NULL) { 
        pred_time = flash_calculation_get_time(NULL) - pred_time;
        split_pred_time += pred_time;

        return 0; 
    }

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

    pred_time = flash_calculation_get_time(NULL) - pred_time;
    split_pred_time += pred_time;

    return 1;
}

int flash_calculation_stab_ann_predict(FLASH_STAB_ANN *fsa, double *input,
        int n, int *stable)
{
    int i;
    FLASH_TENSOR *ft;
    double Pu, Pd, P;
    double pred_time;

    pred_time = flash_calculation_get_time(NULL);

    if (fsa->ann_upper == NULL) {
        pred_time = flash_calculation_get_time(NULL) - pred_time;
        stab_pred_time += pred_time;

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

    pred_time = flash_calculation_get_time(NULL) - pred_time;
    stab_pred_time += pred_time;

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
    FLASH_STAB_ANN *fsa0 = *fsa;

    if (fsa0 != NULL && fsa0->ann_upper != NULL) {
        flash_calculation_ann_model_free(&(fsa0->ann_upper));
        flash_calculation_ann_model_free(&(fsa0->ann_down));
    }

    free(fsa0->ann_upper);
    free(fsa0->ann_down);

    free(*fsa);
}

double flash_calculation_split_pred_time_cost()
{
    double split_pred_time0 = split_pred_time;

    MPI_Allreduce(&split_pred_time0, &split_pred_time, 
            1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return split_pred_time;
}


double flash_calculation_stability_pre_time_cost()
{
    double stab_pred_time0 = stab_pred_time;

    MPI_Allreduce(&stab_pred_time0, &stab_pred_time, 
            1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return stab_pred_time;
}

