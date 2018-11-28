DEVICEFUNC
void polarization_vector(double k[4], sim5complex wp, sim5metric *metric, double f[4])
//! Photon polarization vector.
//! The returned polarization vector satisfies f.k=0 and f.f=1. Since f can be freely shifted 
//! by a multiple of k (f' = f + \alpha*k), it has a freedom in one compoment and it can always 
//! be set in such a way that its time-component is zero (f[0]=0). This routine returns f in such a form.
//!
//! @param k photon 4-momentum vector
//! @param wp complex Walker-Penrose constant
//! @param m metric
//! @param f photon polarization vector (output)
//!
//! @result Photon polarization vector `f`.
{
    /*
    // the followin is the Dexter2016 version, which does not seem to work (k.g > 0)
    double a = metric->a;
    double m = metric->m;
    double r = metric->r;
    double r2 = r*r;
    double a2 = a*a;
    double s2 = 1.0-m*m;
    double s = sqrt(s2);

    // Dexter(2016), Eq 17-22
    double d1 = r*k[0] - r*a*s2*k[3];
    double d2 = a2*s*m*k[0] - a*s*m*(r2+a2)*k[3];
    double d3 = a*r*s2*k[1] + a*s*m*(r2+a2)*k[2];
    double g1 = a*m*k[0] - a2*m*s2*k[3];
    double g2 = -a*r*s*k[0] + r*(r2+a2)*s*k[3];
    double g3 = a2*m*s2*k[1] - r*(r2+a2)*s*k[2];
    double K1 = creal(wp);
    double K2 = cimag(wp);
    
    // Dexter(2016), Eq 23-28
    f[0] = 0.0;
    f[1] =  (g2*K1-d2*K2)*(metric->g33*k[3]+metric->g03*k[0]) - (g3*K1-d3*K2)*metric->g22*k[2];
    f[2] = -(g1*K1-d1*K2)*(metric->g33*k[3]+metric->g03*k[0]) - (g3*K1-d3*K2)*metric->g11*k[1];
    f[3] =  (g1*K1-d1*K2)*(metric->g22*k[2])                  - (g2*K1-d2*K2)*metric->g11*k[1];

    vector_norm_to(f, 1.0, metric);
    */

    // the followin is the original version, which has been replaced by the code above to keep the consistency
    // with Dexter's definitions. Perhaps the two codes give the same result, but that should be tested
    double a = metric->a;
    double m = metric->m;
    double r = metric->r;
    double s = sqrt(1.0-m*m);
    double ra2 = r*r + a*a;
    double r2 = r*r;
    double a2 = a*a;
    double s2 = 1.0-m*m;

    // assert that theta>0.0
    if (s < 1e-12) {
        s  = 1e-12;
        s2 = 1e-24;
        m  = 1.0-0.5*s;
    }

    // obtain A1, A2 from the definition of Walker-Penrose constant
    // see photon_wp_const() function
    double A1 = (+r*creal(wp) - a*m*cimag(wp))/(r*r + a*a*m*m);
    double A2 = (-r*cimag(wp) - a*m*creal(wp))/(r*r + a*a*m*m);

    // f can be freely shifted by a multiple of k (f' = f + \alpha*k), therefore
    // it has a freedom in one compoment and we can always set f in such a way
    // that its time-component will be zero (f[0]=0)
    // then we can solve this set of equations
    //   A1 = (k[0]*f[1]-k[1]*f[0]) + a*s*s*(k[1]*f[3]-k[3]*f[1]);
    //   A2 = s*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
    // along with the orthogonality condition
    //   k.f = 0
    // and we will obtain (Maple was used for that) the following for components of f
    // note: using k.f=0 gives much simpler result than using f.f=1 for the third equation instead
    f[0] = 0.0;
    f[3] = (
              + metric->g11*A1*k[1]*(s*r2*k[3] + s*a2*k[3] - s*a*k[0])
              + metric->g22*A2*k[2]*(k[0] - a*s2*k[3])
           ) / (
             + sqr(k[0])*metric->g33*(s*k[3]*a)
             + sqr(k[0])*metric->g03*(s*k[0]*a - s*r2*k[3] - s*a2*k[3] - a2*s*s2*k[3])
             + sqr(k[1])*metric->g11*a*s*s2*(+ r2*k[3] + a2*k[3] - a*k[0])
             + sqr(k[2])*metric->g22*(a2*a*s*s2*k[3] + r2*a*s*s2*k[3] - s*r2*k[0] - s*a2*k[0])
             + sqr(k[3])*metric->g33*s*(k[3]*a*s2*r2 + k[3]*a2*a*s2 - k[0]*r2 - k[0]*a2 - a2*s2*k[0])
             + sqr(k[3])*metric->g03*a*s*s2*(r2*k[0] + a2*k[0])
           );
    f[1] = (A1-a*s*s*k[1]*f[3])/(k[0]-a*s*s*k[3]);
    f[2] = (A2 + s*k[2]*f[3]*ra2)/(s*k[3]*ra2 - s*a*k[0]);

    vector_norm_to(f, 1.0, metric);
}



/*
DEVICEFUNC
sim5complex photon_wp_const(double k[4], double f[4], sim5metric *metric)
//! Walker-Penrose constant of null geodesic.
//! Calculates components of Walker-Penrose (Walker&Penrose 1970) constant following 
//! Connors, Piran, Stark (1980): 
//! kappa_wp = kappa1 + I*kappa2 = (A1 - I*A2)*(r - I*a*cos(theta)) => 
//! kappa1 = +r*A1 - a*cos(theta)*A2; kappa2 = -r*A2 - a*cos(theta)*A1
//! returns kappa_wp = kappa1 + I*kappa2
//! Note the definition of kappa1 & kappa2, which is opposite to CPS(1980) and 
//! is the same as Dexter(2016) Eq.(3).
//!
//! @param k photon 4-momentum vector (satisfies k.k=0)
//! @param f photon polarization vector (satisfies f.k=0)
//! @param m metric
//!
//! @result Complex Walker-Penrose constant.
{
    double a = metric->a;
    double r = metric->r;
    double m = metric->m;
    double s = sqrt(1.0-m*m);

    // evaluate Walker-Penrose constant components
    // this is the same as Dexter(2017) Eq.(3)
    double A1 = (k[0]*f[1]-k[1]*f[0]) + a*s*s*(k[1]*f[3]-k[3]*f[1]);
    double A2 = s*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
    double wp1 = +r*A1 - a*m*A2;
    double wp2 = -r*A2 - a*m*A1;

    return wp1 + ComplexI*wp2;
}
*/


DEVICEFUNC
sim5complex polarization_constant(double k[4], double f[4], sim5metric *metric)
//! Walker-Penrose constant for a null geodesic.
//! Calculates components of Walker-Penrose (Walker&Penrose 1970) constant 
//! following Dexter(1916).
//! returns 
//! Note: K1 and K2 relate to Connors,Piran,Stark(1980) kappa1=K2, kappa2=K1.
//!
//! @param k photon 4-momentum vector (satisfies k.k=0)
//! @param f photon polarization vector (satisfies f.k=0)
//! @param metric local metric
//!
//! @result Complex Walker-Penrose constant $K_{wp} = K1 + I*K2$.
{
    double a = metric->a;
    double m = metric->m;
    double r = metric->r;

    // following Connors, Piran, Stark (1980):
    double A1 = (k[0]*f[1]-k[1]*f[0]) + a*(1.-m*m)*(k[1]*f[3]-k[3]*f[1]);
    double A2 = sqrt(1.-m*m)*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
    double wp1 = +r*A1 - a*m*A2;
    double wp2 = -r*A2 - a*m*A1;

    return wp1 + ComplexI*wp2;

/*
    double a = metric->a;
    double m = metric->m;
    double r = metric->r;
    double s = sqrt(1.0-m*m);
    double r2 = r*r;
    double a2 = a*a;
    double s2 = 1.0-m*m;

    // Dexter(2016), Eq 17-22
    double d1 = r*k[0] - r*a*s2*k[3];
    double d2 = a2*s*m*k[0] - a*s*m*(r2+a2)*k[3];
    double d3 = a*r*s2*k[1] + a*s*m*(r2+a2)*k[2];
    double g1 = a*m*k[0] - a2*m*s2*k[3];
    double g2 = -a*r*s*k[0] + r*(r2+a2)*s*k[3] ;
    double g3 = a2*m*s2*k[1] - r*(r2+a2)*s*k[2];
    
    // Dexter(2016), Eq 14-15
    double K1 = d1*f[1] + d2*f[2] + d3*f[3];
    double K2 = g1*f[1] + g2*f[2] + g3*f[3];

    return K1 + ComplexI*K2;
*/


    /*
    // after Connors, Piran, Stark (1980): 
    // [same as Dexter(2017) Eq.(3)]
    // K_wp = K1 * i*K2 = (A1 - I*A2)*(r - I*a*cos(theta)) => 
    // K1 = +r*A1 - a*cos(theta)*A2; 
    // K2 = -r*A2 - a*cos(theta)*A1
    //TODO: check if bellow is the same as above (should be), then bellow is more efficient
    double a = metric->a;
    double r = metric->r;
    double m = metric->m;
    double s = sqrt(1.0-m*m);

    // evaluate Walker-Penrose constant components
    double A1 = (k[0]*f[1]-k[1]*f[0]) + a*s*s*(k[1]*f[3]-k[3]*f[1]);
    double A2 = s*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
    double wp1 = +r*A1 - a*m*A2;
    double wp2 = -r*A2 - a*m*A1;

    return wp1 + ComplexI*wp2;
    */
}



/*
DEVICEFUNC
sim5complex photon_wp_const_infinity(double a, double inc, double alpha, double beta)
//! Walker-Penrose constant of null geodesic at infinity.
//! Calculates components of Walker-Penrose constant following Dexter(2016), Eq.3 
//! in the asymptotic limit of k at infinity - see end of section 2.2.0 in Dexter(2016).
//! kappa_wp = kappa1 + I*kappa2
//! The orientation if the polarization basis (at infinity) is that Q=1 corresponds to 
//! horizontal polarization, i.e. the polarizaion vector is along local x-axis in the plane 
//! orthogonal to k (along the the local phi-vector).
//!
//! @param a        black-hole spin
//! @param inc      observer's inclination [rad]
//! @param alpha    alpha impact paramatter
//! @param beta     beta impact parameter
//!
//! @result Complex Walker-Penrose constant.
{
    // X vector (e_phi):   wp1=-gamma, wp2=-beta
    // Y vector (e_theta): wp1=-beta,  wp2=+gamma
    // gamma = -alpha - a*sim(theta0)

    // WP constant for X vector
    double wp1 = -(-alpha - a*sin(inc));
    double wp2 = -beta;
    return wp1 + ComplexI*wp2;
}
*/

DEVICEFUNC
sim5complex polarization_constant_infinity(double a, double alpha, double beta, double incl)
//! Walker-Penrose constant for a null geodesic.
//! Calculates components of Walker-Penrose (Walker&Penrose 1970) constant 
//! following Dexter(1916).
//! returns 
//! Note: K1 and K2 relate to Connors,Piran,Stark(1980) kappa1=K2, kappa2=K1.
//!
//! @param k photon 4-momentum vector (satisfies k.k=0)
//! @param f photon polarization vector (satisfies f.k=0)
//! @param metric local metric
//!
//! @result Complex Walker-Penrose constant $K_{wp} = K1 + I*K2$.
{
    double gamma = - alpha - a*sin(incl);

    double K1 = -gamma;
    double K2 = -beta;

    return K1 + ComplexI*K2;
}


DEVICEFUNC
double polarization_angle_rotation(double a, double inc, double alpha, double beta, sim5complex kappa)
{
// following Connors, Piran, Stark (1980): (note the opposite definition of kappa1,kappa2)
//    double S = l/sin(inc) - a*sin(inc) = -alpha - a*sin(inc)
//    double T = sqrt(q - pow(l*cos(inc)/sin(inc),2.) + pow(a*cos(inc),2.)) = +/- beta
//    beta>0 for k_\theta<0, beta<0 for k_\theta>0
    double kappa1 = creal(kappa);
    double kappa2 = cimag(kappa);
    double S = -alpha - a*sin(inc);
    double T = +beta;
    double X = (-S*kappa2 - T*kappa1)/(S*S+T*T);
    double Y = (-S*kappa1 + T*kappa2)/(S*S+T*T);
    return atan2(Y,X);
}




/*
DEVICEFUNC
void kappa_pw(double a, double r, double m, double k[4], double f[4], double *kappa1, double *kappa2)
// following Connors, Piran, Stark (1980): kappa_pw = kappa1 - I*kappa2 = (A1 - I*A2)*(r - I*a*cos(theta))
// => kappa1 = +r*A1 - a*cos(theta)*A2; kappa2 = -r*A2 - a*cos(theta)*A1
// returns kappa_pw = kappa1 + I*kappa2
// !! note the definition of kappa1 & kappa2, which is opposite to CPS(1980)
{
    // following Connors, Piran, Stark (1980):
    double A1 = (k[0]*f[1]-k[1]*f[0]) + a*(1.-m*m)*(k[1]*f[3]-k[3]*f[1]);
    double A2 = sqrt(1.-m*m)*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
    *kappa1 = +r*A1 - a*m*A2;
    *kappa2 = -r*A2 - a*m*A1;
}


DEVICEFUNC
double polarization_angle_infty(double a, double inc, double alpha, double beta, sim5complex kappa)
{
// following Connors, Piran, Stark (1980): (note the opposite definition of kappa1,kappa2)
//    double S = l/sin(inc) - a*sin(inc) = -alpha - a*sin(inc)
//    double T = sqrt(q - pow(l*cos(inc)/sin(inc),2.) + pow(a*cos(inc),2.)) = +/- beta
//    beta>0 for k_\theta<0, beta<0 for k_\theta>0
    double kappa1 = creal(kappa);
    double kappa2 = cimag(kappa);
    double S = -alpha - a*sin(inc);
    double T = +beta;
    double X = (+S*kappa2 - T*kappa1)/(S*S+T*T);
    double Y = (-S*kappa1 + T*kappa2)/(S*S+T*T);
    return atan2(Y,X);
}
*/


/*

DEVICEFUNC INLINE
void stokes_add(stokesp* dest, stokesp value)
{
    dest->i += value.i;
    dest->q += value.q;
    dest->u += value.u;
    dest->v += value.v;
}


DEVICEFUNC INLINE
double stokes_poldeg(stokesp sp)
{
    return (sp.i>0.0) ? sqrt(sqr(sp.q) + sqr(sp.u))/sp.i : 0.0;
}


DEVICEFUNC INLINE
double stokes_polang(stokesp sp)
{
    double pxang = (sp.i>0.0) ? 0.5*atan2(sp.u/sp.i,sp.q/sp.i) : 0.0;
    while (pxang < 0.0)  pxang += 2.*M_PI;
    while (pxang > M_PI) pxang -= M_PI;
    return pxang;
}


DEVICEFUNC INLINE
void lorentz_boost_x(double V, double X[])
// makes lorentz boost of vector X(t,x,y,z) in x direction by velocity V
{
    double t = X[0];
    double x = X[1];
    double gamma = 1./sqrt(1.-V*V);
    X[0] = gamma * (t + V*x);
    X[1] = gamma * (x + V*t);
}

*/

