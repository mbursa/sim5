//************************************************************************
//    sim5kerr.h - data file
//------------------------------------------------------------------------
//    Date   : 12.10.2004
//    Author : Michal Bursa
//    e-mail : bursa@sirrah.troja.mff.cuni.cz
//------------------------------------------------------------------------
//    (C) 2004 Michal Bursa
//************************************************************************
// $Id: sim4datafile.h,v 1.7 2005/12/08 10:36:12 bursa Exp $

#ifndef _SIM5KERR_GEOD_H
#define _SIM5KERR_GEOD_H

typedef struct geodesic {
	double alpha;                            // impact alpha (x)
	double beta;                             // impact beta (y)
	double i;                                // inclination
	double cos_i;                            // cosine of inclination
	double a;                                // spin
	double l;                                // L_z / E_\infty
	double q;                                // L / E_\infty^2
	complex r1,r2,r3,r4;                     // roots of R-integral
	int    nrr;                              // number of real roots of R integral
	double m2p,m2m,m2,mK;                    // roots and coeficients of T-integral
	double rp;                               // radial distance of periastron, or 0.0 if no periastron exists (no real roots)
	double rp_R;                             // value of R-integral at periastron, or 0.0 if no periastron exists (no real roots)
	double k[4];                             // photon 4-momentum
	double x;                                // value of R/T integrals for some particular solution (code's affine parameter)
//	double _dxds;
//	double _t;
} geodesic;


DEVICEFUNC void geodesic_s2i_setup(double i, double a, double alpha, double beta, geodesic *g, int *status);
DEVICEFUNC void geodesic_s2i_setup_local(double a, double l, double q, geodesic *g, int *status);

DEVICEFUNC double geodesic_s2i_int_R_periastron(geodesic *g);
DEVICEFUNC double geodesic_s2i_int_R(geodesic *g, double r, int beyond_pa);
DEVICEFUNC double geodesic_s2i_inv_R(geodesic *g, double R, int* beyond_pa);

DEVICEFUNC double geodesic_s2i_int_T_eqplane(geodesic *g, int order, double *dk2dx);
DEVICEFUNC double geodesic_s2i_inv_T(geodesic *g, double T, double *dk2dx);

//DEVICEFUNC double geodesic_s2i_polar_eqplane(geodesic *g, double r, int beyond_pa);

DEVICEFUNC double geodesic_s2i_phi(geodesic *g, double x, double *opt_r, double *opt_m);
DEVICEFUNC double geodesic_s2i_timedelay(geodesic *g, double x, double *opt_r, double *opt_m);
DEVICEFUNC double geodesic_s2i_afp(geodesic *g, double x, double *opt_r, double *opt_m);

DEVICEFUNC void geodesic_s2i_solution_eqplane(geodesic *g, int order, double *r, int *beyond_pa, double *phi, double *time_delay, int *status);
DEVICEFUNC void geodesic_s2i_solution_surface(geodesic *g, double *r, double *m, double (*H)(double), double accuracy, int *status);


//DEVICEFUNC double int_theta(double i, double a, double l, double q, double alpha, double beta);
//DEVICEFUNC double _R(double r, double l, double q, double a);
//DEVICEFUNC double int_r(double i, double a, double l, double q, double alpha, double beta);

DEVICEFUNC void geodesic_s2i_follow_init(geodesic *g, double rmax, int *status);
//DEVICEFUNC void geodesic_s2i_follow(geodesic *g, double step, double *r, double *m, double *phi, double *t, double *ds, int *status);
DEVICEFUNC void geodesic_s2i_follow(geodesic *g, double step, double *r, double *m, double *phi, double *t, int *status);
 

#endif
