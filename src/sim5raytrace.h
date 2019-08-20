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

#ifdef __cplusplus
extern "C" {
#endif


// options for raytracing
#define RTOPT_NONE              0         // option for default actions
#define RTOPT_FLAT              1         // assumes Minkowski (flat) metric instead of Kerr
#define RTOPT_POLARIZATION      2         // track change of polarization vector


typedef struct raytrace_data {
    int opt_gr;             // the metric: 1=Kerr metric, 0=flat metric
    int opt_pol;            // polarization: 1=follow transport of f, 0=ignore f
    double step_epsilon;    // step size control factor (note: precision ~ step^2)

    double bh_spin;         // black hole spin
    double E;               // initial motion constant - energy (k_t)
    double Q;               // initial motion constant - Carter constant
    sim5complex WP;         // Walker-Penrose constant (K_wp = wp1 + i*wp2)

    // runtime variables    
    int pass;               // number of passes to raytrace() routine (steps along geodesics)
    int refines;            // number of refinements (recursive calls within raytrace())
    double dk[4];           // derivative of momentum vector in the last step
    double df[4];           // derivative of polarization vector in the last step
    double kt;              // motion constant k_t in the last step
    float error;            // fractional error in the last step
} raytrace_data;


DEVICEFUNC
void raytrace_prepare(double bh_spin, double x[4], double k[4], double presision_factor, int options, raytrace_data* rtd);

DEVICEFUNC
void raytrace(double x[4], double k[4], double *step, raytrace_data* rtd);

DEVICEFUNC
double raytrace_error(double x[4], double k[4], raytrace_data* rtd);

#ifdef __cplusplus
}
#endif


#endif
