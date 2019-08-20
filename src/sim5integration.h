//************************************************************************
//    sim5integration.h - integration functions
//------------------------------------------------------------------------
//    Written by:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1a, 141-31 Praha 4, Czech Republic
//------------------------------------------------------------------------
//    Created: 09.12.2015
//************************************************************************



#ifndef _SIM5INTEGRATION_H
#define _SIM5INTEGRATION_H

#ifdef __cplusplus
extern "C" {
#endif

DEVICEFUNC double integrate_trapezoid(double(*f)(double), double a, double b, double acc);
DEVICEFUNC double integrate_simpson(double (*f)(double), double a, double b, double acc);

DEVICEFUNC void gauleg(double x1, double x2, double x[], double w[], int n);

#ifdef __cplusplus
}
#endif

#endif
