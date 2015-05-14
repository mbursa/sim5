//************************************************************************
//    SIM5 library
//    sim5raytrace.h - photon geodesic integration (direct method)
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************


#ifndef _SIM5_RAYTRACE_H
#define _SIM5_RAYTRACE_H

typedef struct raytrace_data {
    long pass;              // number of passes

    double Q;               // Carter constant for the initial direction
    double K1,K2;           // Wolker-Penrose constant (Kwp = K1 + iK2)
    
    double kt;              // motion constant k_t in the last step
    double kf;              // motion constant k_\phi in the last step
    float error;            // fractional error in the last step
} raytrace_data;


DEVICEFUNC
void kerr_raytrace_prepare(double bh_spin, double x[4], double k[4], double dk[4], double f[4], raytrace_data* rtd);

DEVICEFUNC
void kerr_raytrace(double bh_spin, double x[4], double k[4], double dk[4], double f[4], double *step, raytrace_data* rtd);


#endif
