//************************************************************************
//    sim5utils.h - data file
//------------------------------------------------------------------------
//    Date   : 12.10.2004
//    Author : Michal Bursa
//    e-mail : bursa@sirrah.troja.mff.cuni.cz
//------------------------------------------------------------------------
//    (C) 2004 Michal Bursa
//************************************************************************
// $Id: sim4datafile.h,v 1.7 2005/12/08 10:36:12 bursa Exp $

#ifndef _SIM5POLYROOTS_H
#define _SIM5POLYROOTS_H


DEVICEFUNC int quadratic_eq(double pr, double pi, double qr, double qi, double zr[2], double zi[2]);
DEVICEFUNC int cubic_eq(double p, double q, double r, double zr[3], double zi[3]);
DEVICEFUNC int quartic_eq(double a3, double a2, double a1, double a0, double zr[4], double zi[4]);
DEVICEFUNC void sort_roots_re(double *r1, double *r2, double *r3, double *r4);
DEVICEFUNC void sort_mix(double *r1, double *r2, int *s);
DEVICEFUNC void sort_mix2(double *r1, double *r2, int *s);
//void sort_roots(int *s, double complex *z1, double complex *z2, double complex *z3, double complex *z4);
DEVICEFUNC void sort_roots(int *s, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4);

//void quartic_eq_c(double a3, double a2, double a1, double a0, int *nr, double complex *z1, double complex *z2, double complex *z3, double complex *z4);
DEVICEFUNC void quartic_eq_c(double a3, double a2, double a1, double a0, int *nr, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4);

#endif
