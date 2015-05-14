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
    sim5metric m;
};
typedef struct sim5tetrad sim5tetrad;



DEVICEFUNC
void flat_metric(double a, double r, double m, sim5metric *metric);

DEVICEFUNC
void kerr_metric(double a, double r, double m, sim5metric *metric);

DEVICEFUNC
void kerr_connection(double a, double r, double m, double G[4][4][4]);

DEVICEFUNC
void flat_connection(double a, double r, double m, double G[4][4][4]);

DEVICEFUNC INLINE
void Gamma(double G[4][4][4], double v[4], double result[4]);

DEVICEFUNC INLINE
void vector_covariant(double V1[4], double V2[4], sim5metric* m);

DEVICEFUNC INLINE
double vector_norm(double V[4], sim5metric* m);

DEVICEFUNC INLINE
void vector_norm_to(double V[4], sim5metric* m, double norm);

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



DEVICEFUNC
void photon_momentum(double a, double r, double m, double l, double q2, double r_sign, double m_sign, double k[4]);

DEVICEFUNC
void photon_motion_constants(double a, double r, double m, double k[4], double* L, double* Q);

DEVICEFUNC
double photon_carter(sim5metric *metric, double k[4]);


