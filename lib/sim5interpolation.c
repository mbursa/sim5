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




DEVICEFUNC INLINE
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
//     if ( x <= x[0] )               then  index == 0
//     if ( x >= x[i] && x < x[i+1] ) then  index == i
//     if ( x >= x[N] )               then  index == N-1
{
    long ilo = index_lo;
    long ihi = index_hi;
    while(ihi > ilo + 1) {
        long i = (ihi + ilo)/2;
        if (x < x_array[i])
            ihi = i;
        else
            ilo = i;
    }
    return ilo;
}




DEVICEFUNC INLINE 
static long sim5_interp_search_accel(sim5interp* interp, double x)
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



DEVICEFUNC
static void spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
//! Calculates second derivatives of function y[]=f(x[]) for cubic spline interpolation.
//! - given arrays x[] and y[] containing a tabulated function y=f(x),
//!   with and given values yp1 and ypn for the first derivative of the interpolating
//!   function at points 1 and n , respectively
//! - the routine returns an array y2[1..n] that contains the second derivatives 
//!   of the interpolating function at the tabulated points x i . 
//! - if yp1 and/or ypn are larger than 1e30, respectively, the routine is signaled 
//!   to set the corresponding boundary condition for a natural spline, with zero 
//!   second derivative on that boundary
//! (routine from Numerical Recipes in C)
{
    int i,k;
    double p, qn, sig, un, *u;

    //MALLOC(u,float,n-1);
    u = (double*)malloc((n-1)*sizeof(double));
    if (u == NULL) exit(EXIT_FAILURE);
    //fprintf(stderr,"ERR %s line %d: Memory allocation failure.\n",  __FILE__, __LINE__);

    if(yp1 > 0.99e30)
        y2[0] = u[0] = 0.0;
    else{
        y2[0] = -0.5;
        u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }

    for(i = 1; i < n-1; i++){
        sig = (x[i] - x[i-1])/(x[i+1] - x[i-1]);
        p = sig*y2[i-1] + 2.0;
        y2[i] = (sig - 1.0)/p;
        u[i] = (y[i+1] - y[i])/(x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1]);
        u[i] = (6.0*u[i]/(x[i+1] - x[i-1]) - sig*u[i-1])/p;
    }

    if(ypn > 0.99e30)
        qn = un = 0.0;
    else{
        qn = 0.5;
        un = (3.0/(x[n-1] - x[n-2]))*(ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]));
    }

    y2[n-1] = (un - qn*u[n-2])/(qn*y2[n-2] + 1.0);

    for(k = n-2; k >= 0; k--){
        y2[k] = y2[k]*y2[k+1] + u[k];
    }

    free(u);
}



DEVICEFUNC
static double splint(double xa[], double ya[], double y2a[], int n, double x)
//! Cubic spline interpolation.
//! - given the arrays xa[] and ya[] of dimension N, which tabulate a function and
//!   given the array y2a[] , which is the output from spline() routine,
//!   this routine returns a cubic-spline interpolated value y at point x
//! - xa[] must be orderd array
//! (routine from Numerical Recipes in C)
{
    int  klo,khi,k;
    double h,b,a;
    static int pklo=0, pkhi=1;

    #pragma omp threadprivate(pklo,pkhi)

    // Assuming that calls to this function are made with closely-spaced, 
    // steadily-increasing values of x, we first try using the same values of klo and khi 
    // as were used in the previous invocation. 
    // If that interval is no longer correct, a standard binary search looks for the correct interval.
    if(xa[pklo] <= x && xa[pkhi] > x){
        klo = pklo;
        khi = pkhi;
    }
    else{
        klo = 0;
        khi = n-1;
        while (khi-klo > 1){
            k = (khi + klo) >> 1;
            if(xa[k] > x) khi = k; else klo = k;
        }
        pklo = klo;
        pkhi = khi;
    }

    h = xa[khi] - xa[klo];
    // we can skip this check since have checked that already during sim5interp initialization
    //if (h == 0.0) {
    //    fprintf(stderr,"-E- %s line %d: Bad xa input to function splint()\n", __FILE__,__LINE__);
    //    exit(EXIT_FAILURE);
    //}
    a = (xa[khi] - x)/h;
    b = (x - xa[klo])/h;
    return  a*ya[klo] + b*ya[khi] + ((a*a*a - a)*y2a[klo] + (b*b*b - b)*y2a[khi])*(h*h)/6.0;
}



DEVICEFUNC
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



DEVICEFUNC
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
    if ((interp_type==INTERP_TYPE_SPLINE) && (interp_options & INTERP_OPT_CAN_EXTRAPOLATE)) {
        fprintf(stderr, "ERR (sim5_interp_init): spline interpolation cannot be used with extrapolation option\n");
        return;
    }    


    interp->datamodel = data_model;
    interp->type      = interp_type;
    interp->options   = interp_options;
    interp->d2Y       = NULL;

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



DEVICEFUNC
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



DEVICEFUNC
double sim5_interp_eval(sim5interp* interp, double x)
// makes the evalutaion if interpolated grid at point x
{
    double x_lo, x_hi;
    double y_lo, y_hi;
    long index;

    // treat spline interpolation seperately    
    if (interp->type == INTERP_TYPE_SPLINE) {
        // calculate second derivatives if they are not yet available
        if (!interp->d2Y) {
            interp->d2Y = (double*) malloc(interp->N*sizeof(double));
            spline(interp->X, interp->Y, interp->N, 1e50, 1e50, interp->d2Y);
        }
        return splint(interp->X, interp->Y, interp->d2Y, interp->N, x);
    }



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

    // seems unnecessary as we have checked on order of X array already on initialization
    //if (x_lo >= x_hi) {
    //    fprintf(stderr, "ERR (sim5_interp_eval: unordered X grid (x[%ld]=%.4e, x[%ld]=%.4e, N=%ld)\n", index, x_lo, index+1, x_hi, interp->N);
    //    return NAN; 
    //}

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


/*
double sim5_interp_integral(sim5interp* interp, double a, double b)
// makes the evalutaion of interpolated grid at point x
{
    int i, N = 500;
    double result = 0.0;
    for (i=0; i<N; i++) result += sim5_interp_eval(interp, a+(i+0.5)*(b-a)/(N));
    return result*(b-a)/N;
}
*/



DEVICEFUNC
void sim5_interp_done(sim5interp* interp)
// function frees the interpolation object interp (including a copied data, if necessary)
{
    if ((interp->datamodel==INTERP_DATA_COPY) || (interp->datamodel==INTERP_DATA_BUILD)){
        free(interp->X);
        free(interp->Y);
    }

    if (interp->d2Y) free(interp->d2Y);

    interp->N = 0;
    interp->capa = 0;
    interp->X = NULL;
    interp->Y = NULL;
}



DEVICEFUNC
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



