//************************************************************************
//    sim5interpolation.c - interpolation functions
//------------------------------------------------------------------------
//    Written by:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1a, 141-31 Praha 4, Czech Republic
//------------------------------------------------------------------------
//    Created: 24.10.2012
//************************************************************************



/*
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "sim5interpolation.h"
*/




inline
long sim5_interp_search(const double x_array[], double x, long index_lo, long index_hi)
// perform a search of an array of values
// - the parameters index_lo and index_hi provide an initial bracket, and it is assumed 
//   that index_lo < index_hi
// - the resulting index is guaranteed to be strictly less than index_hi and greater than
//   or equal to index_lo, so that the implicit bracket [index, index+1] always corresponds
//   to a region within the implicit value range of the value array
//   (note that this depends on the result region, i.e. the behaviour at the boundaries 
//   may not correspond to what you expect)
// - complete specification of the behaviour is the following:
//     suppose the input x_array[] = { x0, x1, ..., xN }
//     if ( x == x0 )           then  index == 0
//     if ( x > x0 && x <= x1 ) then  index == 0, and sim. for other interior pts
//     if ( x == xN )           then  index == N-1
//     if ( x > xN )            then  index == N-1
//     if ( x < x0 )            then  index == 0 
{
    long ilo = index_lo;
    long ihi = index_hi;
    while(ihi > ilo + 1) {
        long i = (ihi + ilo)/2;
        if(x_array[i] > x)
            ihi = i;
        else
            ilo = i;
    }
    return ilo;
}




inline 
long sim5_interp_search_accel(sim5interp* interp, double x)
// performs an accelerated search of an array of values having cached the last used index
{
    long x_index = interp->last_index;
    if(x < interp->X[x_index]) {
        //interp->miss_count++;
        interp->last_index = sim5_interp_search(interp->X, x, 0, x_index);
    } else
    if(x >= interp->X[x_index + 1]) {
        //interp->miss_count++;
        interp->last_index = sim5_interp_search(interp->X, x, x_index, interp->N-1);
    } //else {
        //interp->hit_count++;
    //}
    return interp->last_index;
}




sim5interp* sim5_interp_alloc()
// makes a memory allocation for interpolation object
// - this function is optional and it is for heap-allocated variant of usage:
//     sim5interp* interp;
//     interp = sim5_interp_alloc();
//     sim5_interp_init(interp, ...);
//     sim5_interp_free(interp);
// - stack-allocated variant can also be used like this:     
//     sim5interp interp;
//     sim5_interp_init(&interp, ...);
//     sim5_interp_done(&interp);
{
    return calloc(1, sizeof(sim5interp));
}




void sim5_interp_init(sim5interp* interp, double xa[], double ya[], long N, int data_model, int interp_type, int interp_options)
// initializion of the interpolation object interp for the data (xa,ya) where xa and ya are arrays of size N
// - data_model determines how data should be handled, see description of constants INTERP_DATA_xxx
// - interp_type determines in which way data will be interpolated, see description of constants INTERP_TYPE_xxx
// - interp_options specifies additional options (a combination of options can be used)
// - by default, the interpolation object (interp) does not save the data arrays xa and ya, only saves pointers; 
//   however it make an independent copy of those arrays if INTERP_OPT_COPYDATA option is given (in that case
//   the original arrays can be modified or freed after calling sim5_interp_init)
// - xa data array is always assumed to be strictly ordered, with increasing x values; the behavior for other arrangements undefined
{
    interp->datamodel = data_model;
    interp->type      = interp_type;
    interp->options   = interp_options;

    // check of order
    if ((interp->datamodel==INTERP_DATA_REF) || (interp->datamodel==INTERP_DATA_COPY)) {
        long i;
        for (i=0; i<N-1; i++) {
            if (xa[i] >= xa[i+1]) {
                fprintf(stderr, "ERR (sim5_interp_init): unordered X grid (x[%ld]=%.4e, x[%ld]=%.4e, N=%ld, opt=%d)\n", i, xa[i], i+1, xa[i+1], N, interp_options);
                interp->N = 0;
                interp->X = NULL;
                interp->Y = NULL;
                exit(-1);//return;
            }
        }
    }

    switch (interp->datamodel) {
        case INTERP_DATA_REF:
            // assign the reference
            interp->N = N;
            interp->capa = 0;
            interp->X = xa;
            interp->Y = ya;
            interp->xmin = interp->X[0];
            interp->xmax = interp->X[N-1];
            interp->last_index = (N-1)/2;
            break;

        case INTERP_DATA_COPY:
            // make copy of arrays
            interp->N = N;
            interp->capa = N;
            interp->X = (double*)malloc(N*sizeof(double));
            interp->Y = (double*)malloc(N*sizeof(double));
            memcpy (interp->X, xa, N*sizeof(double));
            memcpy (interp->Y, ya, N*sizeof(double));
            interp->xmin = interp->X[0];
            interp->xmax = interp->X[N-1];
            interp->last_index = (N-1)/2;
            break;

        case INTERP_DATA_BUILD:
            // make copy of arrays
            interp->N = 0;
            interp->capa = N>0 ? N : 100;
            interp->X = (double*)calloc(interp->capa, sizeof(double));
            interp->Y = (double*)calloc(interp->capa, sizeof(double));
            interp->xmin = 0.0;
            interp->xmax = 0.0;
            interp->last_index = 0;
            break;

        default:
            fprintf(stderr, "ERR (sim5_interp_init): unimplemented data model (%d)\n", interp->datamodel);
            exit(-1);//return;
    }

}



void sim5_interp_data_push(sim5interp* interp, double x, double y)
// pushed the pair [x,y] into interpolation array
// (x-data must be orderd)
{
    if (interp->datamodel != INTERP_DATA_BUILD) {
        fprintf(stderr, "ERR (sim5_interp_data_push): you can only push in data with INTERP_DATA_BUILD data model\n");
        return;
    }

    long i = interp->N;

    if ((i>0) && (interp->X[i-1] >= x)) {
        fprintf(stderr, "ERR (sim5_interp_data_push): unordered X grid (x[%ld]=%.4e, x[%ld]=%.4e)\n", i-1, interp->X[i-1], i, x);
        exit(-1);//return;
    }

    interp->X[i] = x;
    interp->Y[i] = y;
    interp->N++;

    if (interp->N >= interp->capa) {
        interp->capa *= 2;
        interp->X = (double*)realloc(interp->X, interp->capa*sizeof(double));
        interp->Y = (double*)realloc(interp->Y, interp->capa*sizeof(double));
    }

    interp->xmin = interp->X[0];
    interp->xmax = interp->X[i];
    interp->last_index = i/2;
}



double sim5_interp_eval(sim5interp* interp, double x)
// makes the evalutaion if interpolated grid at point x
{
    double x_lo, x_hi;
    double y_lo, y_hi;
    long index;

    if ((!(interp->options & INTERP_OPT_CAN_EXTRAPOLATE)) && ((x < interp->xmin) || (x > interp->xmax))) {
        fprintf(stderr, "WRN (sim5_interp_eval): unwarranted extrapolation (x=%.4e, xmin=%.4e, xmax=%.4e)\n", x, interp->xmin, interp->xmax);
    }

    if (interp->options & INTERP_OPT_ACCEL) {
        // index search with acceleration
        index = sim5_interp_search_accel(interp, x);
    } else {
        // index search without acceleration
        index = sim5_interp_search(interp->X, x, 0, interp->N-1);
    }


        x_lo = interp->X[index];
        x_hi = interp->X[index + 1];
        y_lo = interp->Y[index];
        y_hi = interp->Y[index + 1];

    if (x_lo >= x_hi) {
        fprintf(stderr, "ERR (sim5_interp_eval: unordered X grid (x[%ld]=%.4e, x[%ld]=%.4e, N=%ld)\n", index, x_lo, index+1, x_hi, interp->N);
        return NAN; 
    }

    switch (interp->type) {
        case INTERP_TYPE_LINLIN:
            return y_lo + (x-x_lo)/(x_hi-x_lo) * (y_hi-y_lo);

        case INTERP_TYPE_LOGLOG:
            return exp(log(y_lo) + (log(x)-log(x_lo)) / (log(x_hi) - log(x_lo)) * (log(y_hi)-log(y_lo)));
            // equvivalent to: exp(log(y_lo) + (log(x)-log(x_lo)) / (log(x_hi) - log(x_lo)) * (log(y_hi)-log(y_lo)))

        case INTERP_TYPE_LOGLIN:
            return y_lo + log(x/x_lo)/log(x_hi/x_lo) * (y_hi-y_lo);
            // equvivalent to: y_lo + (log(x)-log(x_lo)) / (log(x_hi) - log(x_lo)) * (y_hi - y_lo)

        default:
            fprintf(stderr, "ERR (sim5_interp_eval): unimplemented interpolation type (%d)\n", interp->type);
            return NAN;
    }
}


double sim5_interp_integral(sim5interp* interp, double a, double b)
// makes the evalutaion if interpolated grid at point x
{
    int i, N = 500;
    double result = 0.0;
    for (i=0; i<N; i++) result += sim5_interp_eval(interp, a+(i+0.5)*(b-a)/(N));
    return result*(b-a)/N;
}



void sim5_interp_done(sim5interp* interp)
// function frees the interpolation object interp (including a copied data, if necessary)
{
    if ((interp->datamodel==INTERP_DATA_COPY) || (interp->datamodel==INTERP_DATA_BUILD)){
        free(interp->X);
        free(interp->Y);
    }

    interp->N = 0;
    interp->capa = 0;
    interp->X = NULL;
    interp->Y = NULL;
}




void sim5_interp_free(sim5interp* interp)
// function frees the interpolation object interp that had been previously alocated by sim5_interp_alloc
{
    sim5_interp_done(interp);
    free(interp);
}





//#define SIM5FILEIO_TESTIG
#ifdef SIM5FILEIO_TESTIG

int main() {
    double X[5] = {1.,2.,3.,4.,5.};
    double Y[5] = {2.,4.,6.,8.,10.};

    sim5interp interp;
    sim5_interp_init(&interp, X, Y, 5, INTERP_TYPE_LINLIN, INTERP_OPT_ALLOW_EXTRAPOLATION+INTERP_OPT_ACCEL);

    double x;
    for (x=0.; x<10.; x+=0.1) printf("%e %e\n", x, sim5_interp_eval(&interp, x));
    sim5_interp_free(&interp);
    return 0;
}

#endif


