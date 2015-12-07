//************************************************************************
//    sim5kerr-geod - geodesic routines
//------------------------------------------------------------------------
//    Date   : 12.10.2004
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (c) Michal Bursa
//************************************************************************


#ifndef _SIM5KERR_GEOD_H
#define _SIM5KERR_GEOD_H

// structure that holds information about geodesics
typedef struct geodesic {
    // parameters for BH and observer
	double a;                                // BH spin
	double alpha;                            // impact alpha (horizontal)
	double beta;                             // impact beta (vertical)
	double i;                                // inclination
	double cos_i;                            // cosine of inclination

    // geodecis parameters
	double l;                                // motion constant: L_z / E_\infty 
	double q;                                // motion constant: L / E_\infty^2 (Carter's constant)
	sim5complex r1,r2,r3,r4;                 // roots of R-integral
	int    nrr;                              // number of real roots of R integral
	double m2p,m2m,m2,mK;                    // roots and coeficients of T-integral
	double rp;                               // radius of radial turning point (periastron), or 0.0 if no turning point exists (no real roots)

	// evaluation of specific motion integrals
	double Rpa;                              // value of R-integral \int[r_pa..infinity] at radial turning point, or 0.0 trajecotry has no turning point (no real roots)
	double Tpp;                              // value of T-integral \int[-\mu_plus..\mu_plus]
	double Tip;                              // value of T-integral \int[cos_i..\mu_plus]

	double k[4];                             // photon 4-momentum
	double x;                                // value of R/T integrals for some particular solution (code's affine parameter)
} geodesic;


DEVICEFUNC int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *status);
DEVICEFUNC int geodesic_init_src(double a, double r, double m, double k[4], int bpa, geodesic *g, int *status);
DEVICEFUNC double geodesic_P_int(geodesic *g, double r, int bpa);
DEVICEFUNC void geodesic_position(geodesic *g, double P, double x[4]);
DEVICEFUNC double geodesic_position_rad(geodesic *g, double P);
DEVICEFUNC double geodesic_position_pol(geodesic *g, double P);
DEVICEFUNC double geodesic_position_azm(geodesic *g, double r, double m, double P);

DEVICEFUNC double geodesic_P_int(geodesic *g, double r, int bpa);

DEVICEFUNC double geodesic_find_midplane_crossing(geodesic *g, int order);

//DEVICEFUNC void geodesic_s2i_setup(double i, double a, double alpha, double beta, geodesic *g, int *status);
//DEVICEFUNC void geodesic_s2i_setup_local (double a, double l, double q, geodesic *g, int *status);
//DEVICEFUNC void geodesic_s2i_setup_local2(double a, double r, double m, double k[4], int bpa, geodesic *g, int *status);

//DEVICEFUNC double geodesic_s2i_theta_infinity(geodesic *g, double r, double m, int bpa);

//DEVICEFUNC double geodesic_s2i_int_R_periastron(geodesic *g);
DEVICEFUNC double geodesic_s2i_int_R(geodesic *g, double r, int beyond_pa);
DEVICEFUNC double geodesic_s2i_inv_R(geodesic *g, double R, int* beyond_pa);

DEVICEFUNC double geodesic_s2i_int_T_eqplane(geodesic *g, int order, double *dk2dx);
DEVICEFUNC double geodesic_s2i_inv_T(geodesic *g, double T, double *dk2dx);
//DEVICEFUNC double geodesic_s2i_inv_T2(geodesic *g, double T);

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
