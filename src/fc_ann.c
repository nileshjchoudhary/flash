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


static double *flash_calculation_tensor_to_array(FLASH_TENSOR *t)
{
    int i, n;
    double *array;

    n = t->nr * t->nc;

    array = malloc(n * sizeof(*array));

    for (i = 0; i < n; i++) {
        array[i] = t->value[i];
    }

    return array;
}

FLASH_ANN * flash_calculation_ann_model_new(char *file_head, int level,
        FLASH_ANN_TRANS trans)
{
    char file[10240];
    int i;
    FLASH_ANN *ann;
    FLASH_TENSOR *f_scale, *t_scale;

    if (level == 0) {
        return NULL;
    }

    ann = malloc(sizeof(*ann));
    ann->level = level;
    ann->W = malloc(level * sizeof(*(ann->W)));
    ann->b = malloc(level * sizeof(*(ann->b)));

    ann->trans = trans;

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

    sprintf(file, "%s-scale-f.csv", file_head);
    f_scale = flash_calculation_tensor_read(file);

    sprintf(file, "%s-scale-t.csv", file_head);
    t_scale = flash_calculation_tensor_read(file);

    ann->feature_scale = flash_calculation_tensor_to_array(f_scale);
    ann->target_scale = flash_calculation_tensor_to_array(t_scale);

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

    free(ann0->feature_scale);
    free(ann0->target_scale);

    free(*ann);
}

static void flash_calculation_scale_feature(FLASH_ANN *ann, 
        FLASH_TENSOR *x)
{
    int i, j;

    for (i = 0; i < x->nr; i++) {
        for (j = 0; j < x->nc; j++) {
            x->value[i * x->nc + j] = (x->value[i * x->nc + j] - ann->feature_scale[2*j])
                / (ann->feature_scale[2*j + 1] - ann->feature_scale[2*j]);
        }
    }
}

static void flash_calculation_scale_back_target(FLASH_ANN *ann,
        FLASH_TENSOR *x)
{
    int i, j;

    if (ann->trans == FLASH_ANN_TRANS_NONE) {
        for (i = 0; i < x->nr; i++) {
            for (j = 0; j < x->nc; j++) {
                x->value[i * x->nc + j] = x->value[i * x->nc + j]
                    * (ann->target_scale[2*j + 1] - ann->target_scale[2*j])
                    + ann->target_scale[2*j];
            }
        }
    }
    else if (ann->trans == FLASH_ANN_TRANS_LOG) {
        for (i = 0; i < x->nr; i++) {
            for (j = 0; j < x->nc; j++) {
                if (j == 0) {
                    x->value[i * x->nc + j] = x->value[i * x->nc + j]
                        * (ann->target_scale[2*j + 1] - ann->target_scale[2*j])
                        + ann->target_scale[2*j];
                }
                else {
                    x->value[i * x->nc + j] = exp(x->value[i * x->nc + j]);
                }
            }
        }
    }
    else if (ann->trans == FLASH_ANN_TRANS_NEW) {
        double tmp;

        for (i = 0; i < x->nr; i++) {
            for (j = 0; j < x->nc; j++) {
                if (j == 0) {
                    x->value[i * x->nc + j] = x->value[i * x->nc + j]
                        * (ann->target_scale[2*j + 1] - ann->target_scale[2*j])
                        + ann->target_scale[2*j];
                }
                else {
                    tmp = x->value[i * x->nc + j] * x->value[i * x->nc + j];

                    tmp = tmp * (log(ann->target_scale[2*j + 1]) - log(ann->target_scale[2*j]))
                        + log(ann->target_scale[2*j]);

                    x->value[i * x->nc + j] = exp(tmp);
                }
            }
        }
    }
}

double * flash_calculation_predict_value_with_ANN(FLASH_ANN *ann, 
        FLASH_TENSOR *input)
{
    int level = ann->level, i;
    double *value;
    FLASH_TENSOR **y, *x, **W, **b; 

    W = ann->W;
    b = ann->b;

    y = malloc(level * sizeof(*y));
    x = input;

    flash_calculation_scale_feature(ann, x);

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

    flash_calculation_scale_back_target(ann, y[level - 1]);

    value = malloc(y[level - 1]->nr * y[level - 1]->nc * sizeof(double));
    for (i = 0; i < y[level - 1]->nr * y[level - 1]->nc; i++) {
        value[i] = y[level - 1]->value[i];
    }
    //printf("value: %e\n", value);

    for (i = 0; i < level; i++) {
        flash_calculation_tensor_free(&(y[i]));
    }
    free(y);

    return value;
}

FLASH_SPLIT_ANN * flash_calculation_split_ann_model_new(char *file_head, int level, 
        char *target_trans, int ncomp)
{
    char file[1024];
    FLASH_SPLIT_ANN *fsa;
    FLASH_ANN_TRANS trans;

    printf("    Creating ANN ...\n");
    fsa = malloc(sizeof(*fsa));
    fsa->nK = ncomp;

    if (strcmp(target_trans, "none") == 0) {
        trans = FLASH_ANN_TRANS_NONE;
    }
    else if (strcmp(target_trans, "log") == 0) {
        trans = FLASH_ANN_TRANS_LOG;
    }
    else if (strcmp(target_trans, "new") == 0) {
        trans = FLASH_ANN_TRANS_NEW;
    }
    else {
        trans = FLASH_ANN_TRANS_NONE;
    }

    sprintf(file, "%s", file_head);
    printf("file: %s\n", file);
    fsa->ann = flash_calculation_ann_model_new(file, level, trans);

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

    sprintf(file, "%s", file_head);
    printf("file: %s\n", file);
    fsa->ann = flash_calculation_ann_model_new(file, level, FLASH_ANN_TRANS_NONE);

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
    double *value;

    pred_time = flash_calculation_get_time(NULL);

    if (fsa->ann == NULL) { 
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

    value = flash_calculation_predict_value_with_ANN(fsa->ann, ft);

    *F = value[0];
    for (i = 0; i < fsa->nK; i++) {
        K[i] = value[i + 1];
    }

    free(value);

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
    double Pu, Pd, P, *Ps;
    double pred_time;

    pred_time = flash_calculation_get_time(NULL);

    if (fsa->ann == NULL) {
        pred_time = flash_calculation_get_time(NULL) - pred_time;
        stab_pred_time += pred_time;

        return 0;
    }

    ft = malloc(sizeof(*ft));
    ft->nr = 1;
    ft->nc = n - 1;
    ft->value = malloc((n - 1) * sizeof(*(ft->value)));

    for (i = 0; i < n - 1; i++) {
        ft->value[i] = input[i];
    }
    P = input[n - 1];

    Ps = flash_calculation_predict_value_with_ANN(fsa->ann, ft);

    Pu = Ps[0];
    Pd = Ps[1];

    free(Ps);

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
    FLASH_SPLIT_ANN *fsa0 = *fsa;

    if (fsa0 != NULL && fsa0->ann !=NULL ) {
        flash_calculation_ann_model_free(&(fsa0->ann));
    }

    free(*fsa);
}

void flash_calculation_stab_ann_model_free(FLASH_STAB_ANN **fsa)
{
    FLASH_STAB_ANN *fsa0 = *fsa;

    if (fsa0 != NULL && fsa0->ann != NULL) {
        flash_calculation_ann_model_free(&(fsa0->ann));
    }

    free(fsa0->ann);

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

