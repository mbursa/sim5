//************************************************************************
//    sim5distributions.h
//------------------------------------------------------------------------
//    Date   : 10.12.2015
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2015 Michal Bursa
//************************************************************************
//************************************************************************



#ifndef _SIM5DISTRIBUTIONS_H
#define _SIM5DISTRIBUTIONS_H

#ifndef CUDA

#ifdef __cplusplus
extern "C" {
#endif 

typedef struct sim5distrib {
    double x_min;
    double x_max;
    double norm;
    sim5interp pdf;                     // probability distribution function
    sim5interp cdf;                     // cummulative distribution function
    sim5interp icd;                     // inverse cummulative distribution function
} sim5distrib;



DEVICEFUNC void distrib_init(sim5distrib* d, double(*pdf)(double), double x_min, double x_max, int N);
DEVICEFUNC void distrib_done(sim5distrib* d);
DEVICEFUNC INLINE double distrib_hit(sim5distrib* d);

#ifdef __cplusplus
}
#endif 

#endif //CUDA

#endif
