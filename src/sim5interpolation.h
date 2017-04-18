//************************************************************************
//    sim5interpolation.h - interpolation functions
//------------------------------------------------------------------------
//    Written by:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1a, 141-31 Praha 4, Czech Republic
//------------------------------------------------------------------------
//    Created: 24.10.2012
//************************************************************************



#ifndef _SIM5INTERPOLATION_H
#define _SIM5INTERPOLATION_H

#ifndef CUDA

// interpolation options
#define INTERP_DATA_REF                 0       // X/Y arrays are referenced from their supplied original arrays
#define INTERP_DATA_COPY                1       // X/Y arrays are not referenced, but they are copied
#define INTERP_DATA_BUILD               2       // X/Y arrays are not passed, they build by calls to sim5_interp_data_push

// interpolation options
#define INTERP_OPT_ACCEL                1       // interpolation will use acceleration (cashing of index values)
#define INTERP_OPT_CAN_EXTRAPOLATE      2       // extrapolation is allowed when a value for an of-out-grid point is requested

// interpolation types
#define INTERP_TYPE_LINLIN              0       // linear interpolation in both X and Y
#define INTERP_TYPE_LINLOG              1       // linear interpolation in X, logarithmic in Y
#define INTERP_TYPE_LOGLIN              2       // logarithmic interpolation in X, linear in Y
#define INTERP_TYPE_LOGLOG              3       // logarithmic interpolation in both X and Y
#define INTERP_TYPE_SPLINE              4       // linear cubic spline interpolation 


typedef struct sim5interp {
    long    N;                  // X/Y array dimension
    long    capa;               // X/Y array capacity (in case of data is not referenced but stored)
    double* X;                  // array of grid points
    double* Y;                  // array of values
    double* d2Y;                // array of second derivatives (for spline interpolation only)
    int     datamodel;          // data model (INTERP_DATA_xxx)
    int     type;               // interpolation type (INTERP_TYPE_xxx)
    int     options;            // interpolation options (INTERP_OPT_xxx)
    double  xmin;               // minimal value in X grid
    double  xmax;               // maximal value in X grid

    // accelerator:
    long last_index;            // last found index
} sim5interp;


DEVICEFUNC sim5interp* sim5_interp_alloc();
DEVICEFUNC void sim5_interp_init(sim5interp* interp, double xa[], double ya[], long N, int data_model, int interp_type, int interp_options);
DEVICEFUNC void sim5_interp_data_push(sim5interp* interp, double x, double y);
DEVICEFUNC double sim5_interp_eval(sim5interp* interp, double x);
//DEVICEFUNC double sim5_interp_integral(sim5interp* interp, double a, double b);
DEVICEFUNC void sim5_interp_done(sim5interp* interp);
DEVICEFUNC void sim5_interp_free(sim5interp* interp);

DEVICEFUNC INLINE long sim5_interp_search(const double x_array[], double x, long index_lo, long index_hi);

#endif // CUDA

#endif

