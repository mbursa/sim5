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

#ifdef __cplusplus
extern "C" {
#endif

#define GEOD_TYPE_RR               40       // four real roots; allowed region r > r1 (endpoints at infinity)
#define GEOD_TYPE_RR_DBL           41       // four real roots, one double root(r3=r4); allowed region r > r1 
#define GEOD_TYPE_RR_BH            42       // four real roots; allowed region r3 < r < r2 (endpoints under horizon)
#define GEOD_TYPE_RC                2       // two real & two complex roots; allowed region r > r1
#define GEOD_TYPE_CC                0       // four complex roots


#define GD_OK                           0
#define GD_ERROR_Q_ZERO                 1
#define GD_ERROR_BOUND_GEODESIC         2
#define GD_ERROR_UNKNOWN_SOLUTION       3
#define GD_ERROR_TYPE_RR_DOUBLE         4
#define GD_ERROR_TYPE_CC                5
#define GD_ERROR_Q_RANGE                7
#define GD_ERROR_MUPLUS_RANGE           8
#define GD_ERROR_MU0_RANGE              9
#define GD_ERROR_MM_RANGE              10
#define GD_ERROR_INCL_RANGE            11
#define GD_ERROR_SPIN_RANGE            12



// structure that holds information about geodesics
typedef struct geodesic {
    // parameters for BH and observer
	double a;                                // BH spin
	double alpha;                            // impact alpha (horizontal)
	double beta;                             // impact beta (vertical)
	double incl;                            // cosine of inclination
	double cos_i;                            // cosine of inclination

    // geodecis parameters
	double l;                                // motion constant: L_z / E_\infty 
	double q;                                // motion constant: L / E_\infty^2 (Carter's constant)
	sim5complex r1,r2,r3,r4;                 // roots of R-integral
	int    nrr;                              // number of real roots of R integral
	int    type;                             // geodesic type (see GEOD_TYPE_XXX constants)
	double m2p,m2m,mm,mK;                    // roots and coeficients of T-integral

	double rp;                               // radius of radial turning point (pericenter); this is either r1 or r2 (for GEOD_TYPE_RR_BH)
	double dmdp_inf;                         // sign of the derivative d(m)/d(P) at infinity (GEOD_TYPE_RR,GEOD_TYPE_RC); this may be-1.0 or +1.0

	// evaluation of specific motion integrals
	double Rpc;                              // value of R-integral \int[r1..infty] (GEOD_TYPE_RR,GEOD_TYPE_RC) or \int[r3..r2] (for GEOD_TYPE_RR_BH)
	double Tpp;                              // value of T-integral \int[-\mu_plus..\mu_plus]
	double Tip;                              // value of T-integral \int[cos_i..\mu_plus]

	double k[4];                             // photon 4-momentum
	double p;                                // value of R/T integrals for some particular solution (code's affine parameter)
} geodesic;


DEVICEFUNC int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *error);
DEVICEFUNC int geodesic_init_src(double a, double r, double m, double k[4], int bpa, geodesic *g, int *error);
DEVICEFUNC double geodesic_P_int(geodesic *g, double r, int bpa);
DEVICEFUNC void geodesic_position(geodesic *g, double P, double x[4]);
DEVICEFUNC double geodesic_position_rad(geodesic *g, double P);
DEVICEFUNC double geodesic_position_pol(geodesic *g, double P);
DEVICEFUNC double geodesic_position_pol_sign_k_theta(geodesic *g, double P);
DEVICEFUNC double geodesic_position_azm(geodesic *g, double r, double m, double P);
DEVICEFUNC double geodesic_P_int(geodesic *g, double r, int bpa);
DEVICEFUNC double geodesic_dm_sign(geodesic *g, double P);
DEVICEFUNC void geodesic_momentum(geodesic *g, double P, double r, double m, double k[]);
DEVICEFUNC double geodesic_find_midplane_crossing(geodesic *g, int order);
DEVICEFUNC void geodesic_follow(geodesic *g, double step, double *P, double *r, double *m, int *status);
DEVICEFUNC double geodesic_timedelay(geodesic *g, double P1, double r1, double m1, double P2, double r2, double m2);

#ifdef __cplusplus
}
#endif

#endif
