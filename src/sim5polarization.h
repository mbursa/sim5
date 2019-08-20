//************************************************************************
//    SIM5 library
//    sim5polarization.h - polarization routines
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************


#ifndef _SIM5_POLARIZATION_H
#define _SIM5_POLARIZATION_H

#ifdef __cplusplus
extern "C" {
#endif

DEVICEFUNC void polarization_vector(double k[4], sim5complex wp, sim5metric *metric, double f[4]);
DEVICEFUNC sim5complex polarization_constant(double k[4], double f[4], sim5metric *metric);
DEVICEFUNC sim5complex polarization_constant_infinity(double a, double alpha, double beta, double incl);

DEVICEFUNC double polarization_angle_rotation(double a, double inc, double alpha, double beta, sim5complex kappa);

#ifdef __cplusplus
}
#endif


#endif
