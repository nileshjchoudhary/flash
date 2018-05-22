#include "fc.h"

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
