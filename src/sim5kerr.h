//************************************************************************
//    SIM5 library
//    sim5kerr.h - basic properties of Kerr spacetime
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************



struct sim5metric {
    double a,r,m;
    double g00;
    double g11;
    double g22;
    double g33;
    double g03;
};
typedef struct sim5metric sim5metric;

struct sim5tetrad {
    double e[4][4];
    sim5metric metric;
};
typedef struct sim5tetrad sim5tetrad;



DEVICEFUNC
void flat_metric(double r, double m, sim5metric *metric);

DEVICEFUNC
void flat_metric_contravariant(double r, double m, sim5metric *metric);

DEVICEFUNC
void kerr_metric(double a, double r, double m, sim5metric *metric);

DEVICEFUNC
void kerr_metric_contravariant(double a, double r, double m, sim5metric *metric);

DEVICEFUNC
void flat_connection(double r, double m, double G[4][4][4]);

DEVICEFUNC
void kerr_connection(double a, double r, double m, double G[4][4][4]);

DEVICEFUNC INLINE
void Gamma(double G[4][4][4], double U[4], double V[4], double result[4]);

DEVICEFUNC INLINE
void vector(double x[4], double x0, double x1, double x2, double x3);

DEVICEFUNC INLINE
void vector_copy(double src[4], double dst[4]);

DEVICEFUNC INLINE
void vector_covariant(double V1[4], double V2[4], sim5metric* m);

DEVICEFUNC INLINE
double vector_norm(double V[4], sim5metric* m);

DEVICEFUNC INLINE
double vector_3norm(double V[4]);

DEVICEFUNC INLINE
void vector_norm_to(double V[4], double norm, sim5metric* m);

DEVICEFUNC
void vector_norm_to_null(double V[4], double V0, sim5metric* m);

DEVICEFUNC INLINE
void vector_multiply(double V[4], double factor);

DEVICEFUNC INLINE
double dotprod(double V1[4], double V2[4], sim5metric* m);

DEVICEFUNC
void tetrad_zamo(sim5metric *m, sim5tetrad *t);

DEVICEFUNC
void tetrad_radial(sim5metric *m, double v_r, sim5tetrad *t);

DEVICEFUNC
void tetrad_azimuthal(sim5metric *m, double Omega, sim5tetrad *t);

DEVICEFUNC
void tetrad_surface(sim5metric *m, double Omega, double vr, double dhdr, sim5tetrad *t);

DEVICEFUNC
void bl2on(double Vin[4], double Vout[4], sim5tetrad* t);

DEVICEFUNC
void on2bl(double Vin[4], double Vout[4], sim5tetrad* t);


//-----------------------------------------------------------------
// orbital motion
//-----------------------------------------------------------------


DEVICEFUNC INLINE
double r_bh(double a);

DEVICEFUNC INLINE
double r_ms(double a);

DEVICEFUNC INLINE
double r_mb(double a);

DEVICEFUNC INLINE
double r_ph(double a);

DEVICEFUNC INLINE
double OmegaK(double r, double a);

DEVICEFUNC INLINE
double ellK(double r, double a);

DEVICEFUNC INLINE
double omega_r(double r, double a);

DEVICEFUNC INLINE
double omega_z(double r, double a);

DEVICEFUNC INLINE
double Omega_from_ell(double ell, sim5metric *m);

DEVICEFUNC INLINE
double ell_from_Omega(double Omega, sim5metric *m);

DEVICEFUNC INLINE
double gfactorK(double r, double a, double l);

DEVICEFUNC
void photon_momentum(double a, double r, double m, double l, double q2, double r_sign, double m_sign, double k[4]);

DEVICEFUNC
void photon_motion_constants(double a, double r, double m, double k[4], double* L, double* Q);

DEVICEFUNC
double photon_carter_const(double k[4], sim5metric *metric);

DEVICEFUNC
sim5complex photon_wp_const(double k[4], double f[4], sim5metric *metric);

DEVICEFUNC
void polarization_vector(double k[4], sim5complex wp, sim5metric *metric, double f[4]);

DEVICEFUNC
double polarization_angle_infty(double a, double inc, double alpha, double beta, sim5complex kappa);


DEVICEFUNC INLINE
void fourvelocity_zamo(sim5metric *m, double U[4]);


DEVICEFUNC INLINE
void fourvelocity_azimuthal(double Omega, sim5metric *m, double U[4]);

DEVICEFUNC INLINE
void fourvelocity_radial(double vr, sim5metric *m, double U[4]);

DEVICEFUNC INLINE
double fourvelocity_norm(double U1, double U2, double U3, sim5metric *m);


