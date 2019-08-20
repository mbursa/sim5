//************************************************************************
//    SIM5 library
//    sim5kerr.c - basic routines related to Kerr metric
//------------------------------------------------------------------------
//    Author: Michal Bursa (bursa@astro.cas.cz)
//    MIT Licence
//************************************************************************


//! \file sim5kerr.c
//! Basic spacetime routines
//! 
//! Provides basic routines for doing physics in Kerr and Minkowski spacetimes  (metric, connection, tetrads, 
//! vector algebra, orbital motion, photon propagation, etc).

//static long prof_N = 0;
//static double prof_t = 0.0;
//
//void sim5kerr_profile() {
//    fprintf(stderr, "sim5kerr_profile: N=%ld t=%.2e t/N=%.3e\n", prof_N, prof_t, prof_t/(double)(prof_N));
//}



//-----------------------------------------------------------------
// metric, tetrads and vectors
//-----------------------------------------------------------------


DEVICEFUNC
void flat_metric(double r, double m, sim5metric *metric)
//! Flat (Minkowski) spacetime metric.
//! Returns covariant Minkowski metric \f$\eta_\mu\nu\f$ in spherical coordinates (t,r,theta,phi).
//!
//! @param r radial coordinate
//! @param m poloidal coordinate \f$m=\cos\theta\f$
//!
//! @result Metric components are returned in `metric` parameter.
{
    metric->a = 0.0;
    metric->r = r;
    metric->m = m;
    metric->g00 = -1.0;
    metric->g11 = +1.0;
    metric->g22 = +r*r;
    metric->g33 = +r*r*(1.-m*m);
    metric->g03 = 0.0;
}



DEVICEFUNC
void flat_metric_contravariant(double r, double m, sim5metric *metric)
//! Flat (Minkowski) spacetime metric (contravariant).
//! Returns contravariant Minkowski metric \f$\eta^\mu\nu\f$ in spherical coordinates (t,r,theta,phi).
//!
//! @param r radial coordinate
//! @param m poloidal coordinate \f$m=\cos\theta\f$
//!
//! @result Metric components are returned in `metric` parameter.
{
    metric->a = 0.0;
    metric->r = r;
    metric->m = m;
    metric->g00 = -1.0;
    metric->g11 = +1.0;
    metric->g22 = +1./(r*r);
    metric->g33 = +1./(r*r)/(1.-m*m);
    metric->g03 = 0.0;
}



DEVICEFUNC
void kerr_metric(double a, double r, double m, sim5metric *metric)
//*********************************************************
//! Kerr spacetime metric.
//! Returns Kerr metric \f$g_\mu\nu\f$ in spherical coordinates (t,r,theta,phi).
//!
//! @param a black hole spin
//! @param r radial coordinate
//! @param m poloidal coordinate \f$m=\cos\theta\f$
//!
//! @result Metric components are returned in `metric` parameter.
{
    double r2  = sqr(r);
    double a2  = sqr(a);
    double m2  = sqr(m);
    double S   = r2 + a2*m2;
    double s2S = (1.0-m2)/S;
    //double D = r2 - 2.*r + a2;
    //double A = sqr(r2+a2) - a2*D*s2;
    metric->a = a;
    metric->r = r;
    metric->m = m;
    metric->g00 = -1. + 2.0*r/S;
    metric->g11 = S/(r2-2.*r+a2); //S/D;
    metric->g22 = S;
    metric->g33 = (r2*r2 + a2*a2*m2 + a2*r2*(1.+m2) + 2.*r*a2*s2S*S)*s2S;
    metric->g03 = -2.*a*r*s2S;
}



DEVICEFUNC
void kerr_metric_contravariant(double a, double r, double m, sim5metric *metric)
//! Kerr spacetime metric (contravariant).
//! Returns contravariant Kerr metric \f$g^\mu\nu\f$ in spherical coordinates (t,r,theta,phi).
//!
//! @param a black hole spin
//! @param r radial coordinate
//! @param m poloidal coordinate \f$m=\cos\theta\f$
//!
//! @result Metric components are returned in `metric` parameter.
{
    double r2  = sqr(r);
    double a2  = sqr(a);
    double m2  = sqr(m);
    double S   = r2 + a2*m2;
    double SD  = S*(r2 - 2.*r + a2);  // = S*D
    //double D = r2 - 2.*r + a2;
    //double A = sqr(r2+a2) - a2*D*s2;
    metric->a = a;
    metric->r = r;
    metric->m = m;
    metric->g00 = -sqr(r2+a2)/SD + a2*(1.-m2)/S;  // =-A/SD
    metric->g11 = (r2-2.*r+a2)/S; // =D/S
    metric->g22 = 1./S;
    metric->g33 = 1./S/(1.-m2) - a2/SD;
    metric->g03 = -2.*a*r/SD;
}



DEVICEFUNC
void kerr_connection(double a, double r, double m, double G[4][4][4])
//! Christoffel symbol components for Kerr metric (\f$Gamma^\mu_\alpha\beta\f$).
//! Returns a matrix of connection coefficients for Kerr metric.
//!
//! (!) NOTE: Christoffel tensor Gamma^i_(jk) is symmetric in lower two indices (j,k).
//! This function only evaluates components of the tensor, where j<k and
//! multiplies the value of these coeficients by a factor of two.
//! In this way, both summation over j=1..4,k=1..4 and over j=1..4,k=j..4
//! will give the same result; however, care must be taken when summing
//! Gamma^i_(jk) U^j V^k where factor 0.5 has to be used: 
//! `0.5*G[i][j][k]*(u[j]*v[k] + u[k]*v[j])`
//!
//! @param a black hole spin
//! @param r radial coordinate
//! @param m poloidal coordinate \f$m=\cos\theta\f$
//!
//! @result Connection coeficients are returned in `G` parameter.
{
    //prof_N++;
    //clock_t start_t, end_t;
    //start_t = clock();

    double s  = sqrt(1.-m*m);
    double cs = s*m;
    double c2 = m*m;
    double s2 = s*s;
    double cc = c2-s2; //cc
    double CC = 8.*c2*c2-8.*c2+1.;   //cos(4*theta)
    double a2 = a*a;
    double a4 = a2*a2;
    double r2 = r*r;
    double r3 = r2*r;
    double R  = pow(a2 + 2.*r2 + a2*cc, 2.);
    double D  = r2 - 2.*r + a2;
    double S  = r2 + a2*c2;
    double S3 = S*S*S;

    memset(G, 0, 4*4*4*sizeof(double));


    G[0][0][1] = 2.0 * (-2.)*(a2 + r2)*(a2 - 2.*r2 + a2*cc)/(D*R);
    G[0][0][2] = 2.0 * (-8.)*a2*r*cs/R;
    G[0][1][3] = 2.0 * (+2.)*a*s2*(a4 - 3.*a2*r2 - 6.*r2*r2 + a2*cc*(a2 - r2))/(D*R);
    G[0][2][3] = -G[0][0][2]*s2*a;

    G[1][0][0] = -D*(a2*c2-r2)/S3;
    G[1][0][3] = 2.0 * (a*D*(a2*c2-r2)*s2)/S3;
    G[1][1][1] = (r*(a2 - r) + a2*(1. - r)*c2)/(D*S);
    G[1][1][2] = 2.0 * (-a2*cs)/S;
    G[1][2][2] = -r*D/S;
    G[1][3][3] = -D*s2*(2.*c2*a2*r3 + r2*r3 + a4*c2*s2 + c2*c2*a4*r - a2*r2*s2)/S3; //-D*s2*(c2*(2.*a2*r*(a2 + r2) + a4*(1. - r)*s2) + r*(-a4 + r2*r2 + a2*(a2 - r)*s2))/S3;

    G[2][0][0] = -2.*a2*r*cs/S3;
    G[2][0][3] = 2.0 * 2.*a*r*cs*(a2 + r2)/S3;
    G[2][1][1] = a2*cs/S/D;
    G[2][1][2] = 2.0 * r/S;
    G[2][2][2] = -G[2][1][1]*D;
    G[2][3][3] = -2.*cs*(3.*a2*a4 + 10.*a4*r + 11.*a4*r2 + 16.*a2*r3 +
            16.*a2*r2*r2 + 8.*r3*r3 + 4.*a2*(a2 + 2.*r2)*D*cc +
            a4*D*CC)/(16.*S3);

    G[3][0][1] = 2.0 * a*(r2 - a2*c2)/(D*S*S);
    G[3][0][2] = 2.0 * (-8.)*a*r*(m/s)/R;
    G[3][1][3] = 2.0 * (a4 + 3.*a4*r - 12.*a2*r2 + 8.*a2*r3 -
            16.*r2*r2 + 8.*r2*r3 + 4.*a2*r*(2.*r2 -r + a2)*cc -
            a4*(1. - r)*CC)/(2.*D*R);
    G[3][2][3] = 2.0 * ((3.*a4 + 8.*a2*r + 8.*a2*r2 + 8.*r2*r2 +
            4.*a2*(2.*r2 -2.*r + a2)*cc + a4*CC)*(m/s))/(2.*R);

    //end_t = clock();
    //prof_t += (end_t - start_t);
}


DEVICEFUNC
void flat_connection(double r, double m, double G[4][4][4])
//! Christoffel symbol components for Minkowski metric (\f$Gamma^\mu_\alpha\beta\f$).
//! Returns a matrix of connection coefficients for Minkowski metric, i.e. Kerr metric 
//! in the limit of M=0 and a=0.
//!
//! (!) NOTE: Christoffel tensor Gamma^i_(jk) is symmetric in lower two indices (j,k).
//! This function only evaluates components of the tensor, where j<k and
//! multiplies the value of these coeficients by a factor of two.
//! In this way, both summation over j=1..4,k=1..4 and over j=1..4,k=j..4
//! will give the same result; however, care must be taken when summing
//! Gamma^i_(jk) U^j V^k where factor 0.5 has to be used: 
//! `0.5*G[i][j][k]*(u[j]*v[k] + u[k]*v[j])`
//!
//! @param r radial coordinate
//! @param m poloidal coordinate \f$m=\cos\theta\f$
//!
//! @result Connection coeficients are returned in `G` parameter.
{
    double s = sqrt(1.-m*m);
    memset(G, 0, 4*4*4*sizeof(double));

    G[1][2][2] = -r;
    G[1][3][3] = -r*s*s;

    G[2][1][2] = 2.0 * 1./r;
    G[2][3][3] = -m*s;

    G[3][1][3] = 2.0 * 1./r;
    G[3][2][3] = 2.0 * m/s;
}

/*
DEVICEFUNC INLINE
void Gamma(double G[4][4][4], double V[4], double result[4])
// *********************************************************
// returns product of summation -G^i_(jk) V^j V^k - useful for parallel transport
// (note the definition of G tensor components in symmetric indices)
// inputs: Christoffel tensor <G>, vector <V>
// output: a vector
{
    int i,j,k;
    for (i=0;i<4;i++) {
        result[i] = 0.0;
        for (j=0;j<4;j++) for (k=j;k<4;k++) result[i] -= G[i][j][k]*V[j]*V[k];
    }
}
*/


DEVICEFUNC INLINE
void Gamma(double G[4][4][4], double U[4], double V[4], double result[4])
//! Change of vector along a trajectory.
//! Returns product of summation `-G^i_(jk) U^j V^k`. This is useful for calculating 
//! parallel transport of vector `V` along a trajectory specified by tangent vector `U`, 
//! as it gives the derivative \f$dV/d\lambda\f$. (Note the definition of G tensor components 
//! in symmetric indices.)
//!
//! @param G connection coefficients
//! @param U tangent vector
//! @param V transported vector
//!
//! @result Change of transported vector \f$dV/d\lambda\f$.
{
    int i,j,k;
    for (i=0;i<4;i++) {
        result[i] = 0.0;
        for (j=0;j<4;j++) for (k=j;k<4;k++) result[i] -= 0.5*G[i][j][k]*(U[j]*V[k] + U[k]*V[j]);
    }
}


DEVICEFUNC INLINE
void vector_set(double x[4], double x0, double x1, double x2, double x3)
//! Make a 4-vector.
//! Returns a 4-vector that has given components (contravarient form \f$X^\mu\f$).
//!
//! @param x0-x3 components of the vector
//!
//! @result Vector is returned in `x` parameter.
{
    x[0] = x0;
    x[1] = x1;
    x[2] = x2;
    x[3] = x3;
}


DEVICEFUNC INLINE
void vector_copy(double src[4], double dst[4])
//! Copy a 4-vector.
//! Copies a 4-vector `src` into `dst`.
//!
//! @param src the source vector
//! @param dst the target vector
//!
//! @result A copy of src vector in dst.
{
    dst[0] = src[0];
    dst[1] = src[1];
    dst[2] = src[2];
    dst[3] = src[3];
}


DEVICEFUNC INLINE
void vector_covariant(double V1[4], double V2[4], sim5metric* m)
//! Covariant version of a vector.
//! Converts a standard (contravariant) form of a vector to covariant form (\f$X^\mu -> X_\mu\f$).
//! Metric `m` can be NULL in which case flat metric is assumed.
//!
//! @param V1 contravariant vector
//! @param V2 covariant vector (output)
//! @param m metric
//!
//! @result Transformed vector is returned in `V2` parameter.
{
    if (m) {
        V2[0] = V1[0]*m->g00 + V1[3]*m->g03;
        V2[1] = V1[1]*m->g11;
        V2[2] = V1[2]*m->g22;
        V2[3] = V1[3]*m->g33 + V1[0]*m->g03;
    } else {
        V2[0] = -V1[0];
        V2[1] = +V1[1];
        V2[2] = +V1[2];
        V2[3] = +V1[3];
    }
}


DEVICEFUNC INLINE
double vector_norm(double V[4], sim5metric* m)
//! Norm of a vector.
//! Returns norm of a vector V, i.e. sqrt(V*V).
//! Metric `m` can be NULL in which case flat metric is assumed.
//! NOTE: only works for space-like vectors, where V*V>0.
//!
//! @param V vector
//! @param m metric
//!
//! @result Norm of a vector.
{
    return sqrt(dotprod(V, V, m));
}



DEVICEFUNC INLINE
double vector_3norm(double V[4])
//! 3-norm of a vector (in flat spacetime).
//! Returns a 3-norm of a vector V (includes only 'spatial' components) in flat (Minkowski) spacetime, i.e. sqrt(\sum V[i]*V[i]), where i=1..3.
//!
//! @param V vector
//! @param m metric
//!
//! @result 3-norm of a vector.
{
    return sqrt(V[1]*V[1] + V[2]*V[2] + V[3]*V[3]);
}



DEVICEFUNC INLINE
void vector_multiply(double V[4], double factor)
//! Multiply a vector.
//! Multiplies the components of a vector by given factor.
//!
//! @param V vector
//! @param factor multiplication factor
//!
//! @result Modified vector V.
{
    V[0] *= factor;
    V[1] *= factor;
    V[2] *= factor;
    V[3] *= factor;
}



DEVICEFUNC INLINE
void vector_norm_to(double V[4], double norm, sim5metric* m)
//! Normalizes a vector to given size.
//! Changes the vector components such that the norm of the vector becomes `|norm|`, 
//! i.e. \f$V.V = \rm{norm}\f$. The vector `V` must be space-like if norm>0, 
//! and it must be time-like if norm<0 (no checks are performed to enfoce that).
//! Metric `m` can be NULL in which case flat metric is assumed.
//! NOTE: vector cannot be null vector (V*V=0), nor norm can be zero; use vector_norm_to_null() in that case.
//!
//! @param V vector to be normalized
//! @param norm normalization the vector should have 
//! @param m metric
//!
//! @result Changes the components of vector `V`.
{
    double N = dotprod(V, V, m);
    V[0] *= sqrt(norm/N);
    V[1] *= sqrt(norm/N);
    V[2] *= sqrt(norm/N);
    V[3] *= sqrt(norm/N);
}


DEVICEFUNC
void vector_norm_to_null(double V[4], double V0, sim5metric* m)
//! Normalizes a null vector to given size.
//! Changes the vector components such that the time-component is set to V0 and the 
//! remaining components are adjusted such the vector is still a null vector (V.V=0).
//! Metric `m` can be NULL in which case flat metric is assumed.
//! NOTE: the original vector must be null vector (V*V=0) and `V0` cannot be zero (no check is performed on that).
//!
//! @param V vector to be normalized
//! @param V0 the new value of vector time-component
//! @param m metric
//!
//! @result Changes the components of vector `V`.
{
    double a,b,c,alpha;
    if (m) {
        a = V[1]*V[1]*m->g11 + V[2]*V[2]*m->g22 + V[3]*V[3]*m->g33;
        b = V0*V[3]*m->g03;  // this actually is b/2, that is why the formula for roots is altered too
        c = V0*V0*m->g00;
        alpha = fmax(-b/a+sqrt(b*b-a*c)/a, -b/a-sqrt(b*b-a*c)/a);
    } else {
        a = V[1]*V[1] + V[2]*V[2] + V[3]*V[3];
        c = -V0*V0;
        alpha = sqrt(-c/a);
    }
    V[0]  = V0;
    V[1] *= alpha;
    V[2] *= alpha;
    V[3] *= alpha;
}



DEVICEFUNC INLINE
double dotprod(double V1[4], double V2[4], sim5metric* m)
//! Scalar product.
//! Calculates scalar product of two 4-vectors (\f$U^\mu V^\nu g_{\mu\nu}\f$).
//! NOTE: if metric is NULL, flat space metric is used
//!
//! @param V1 first vector
//! @param V1 second vector
//! @param m metric
//!
//! @result The value of scalar product
{
    if (m)
        return V1[0]*V2[0]*m->g00 + V1[1]*V2[1]*m->g11 + V1[2]*V2[2]*m->g22 +
               V1[3]*V2[3]*m->g33 + V1[0]*V2[3]*m->g03 + V1[3]*V2[0]*m->g03;
    else
        return -V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2] + V1[3]*V2[3];
}



DEVICEFUNC
void tetrad_general(sim5metric *m, double U[], sim5tetrad *t)
//! Tetrad of an observer that is moving with a general 4-velocity U.
//! Note: the orientation of theta-vector points in opposite direction than d/d(theta),
//! i.e. at the equatorial plane the theta-vector points in upward direction.
//! (based on Kulkarni+2011, Dexter2016 eq.36-43)
//!
//! @param m metric
//! @param U observer 4-velocity
//! @param t tetrad
//!
//! @result Returns tetrad vectors in `t`.
{
    // shortcuts
    double u[4];
    double D = sqr(m->r) - 2*(m->r) + sqr(m->a);
    vector_covariant(U, u, m);  // U=U^\mu, u = U_\mu

    // norms
    double N1 = sqrt( -m->g11*(u[0]*U[0]+u[3]*U[3])*(1.+u[2]*U[2]) );
    double N2 = sqrt( +m->g22*(1.+u[2]*U[2]) );
    double N3 = sqrt( -(u[0]*U[0]+u[3]*U[3])*D*(1.-sqr(m->m)) );
    
    t->e[0][0] = U[0];
    t->e[0][1] = U[1];
    t->e[0][2] = U[2];
    t->e[0][3] = U[3];

    t->e[1][0] = u[1]*U[0] / N1;
    t->e[1][1] = -(u[0]*U[0]+u[3]*U[3]) / N1;
    t->e[1][2] = 0.0;
    t->e[1][3] = u[1]*U[3] / N1;

    t->e[2][0] = u[2]*U[0] / N2;
    t->e[2][1] = u[2]*U[0] / N1;
    t->e[2][2] = (1.+u[2]*U[2]) / N2;
    t->e[2][3] = u[2]*U[3] / N2;

    t->e[3][0] = -u[0] / N3;
    t->e[3][1] = 0.0;
    t->e[3][2] = 0.0;
    t->e[3][3] = +u[3] / N3;

    t->metric = *m;
}



DEVICEFUNC
void tetrad_zamo(sim5metric *m, sim5tetrad *t)
//! Tetrad of a zero angular momentum observer (ZAMO).
//! Returns basis vectors for ZAMO (locally non-rotating) observer tetrad e_{(a)}^{\mu}.
//! Note the orientation of theta-vector that points in opposite direction than d/d(theta),
//! i.e. at the equatorial plane the theta-vector points in upward direction.
//!
//! @param m metric
//! @param t tetrad
//!
//! @result Returns tetrad vectors in `t`.
{
    t->e[0][0] = sqrt(m->g33/(sqr(m->g03) - m->g33*m->g00));
    t->e[0][1] = 0.0;
    t->e[0][2] = 0.0;
    t->e[0][3] = -t->e[0][0] * m->g03/m->g33;

    t->e[1][0] = 0.0;
    t->e[1][1] = 1./sqrt(m->g11);
    t->e[1][2] = 0.0;
    t->e[1][3] = 0.0;

    t->e[2][0] = 0.0;
    t->e[2][1] = 0.0;
    t->e[2][2] = -1./sqrt(m->g22);
    t->e[2][3] = 0.0;

    t->e[3][0] = 0.0;
    t->e[3][1] = 0.0;
    t->e[3][2] = 0.0;
    t->e[3][3] = 1./sqrt(m->g33);

    t->metric = *m;
}



DEVICEFUNC
void tetrad_radial(sim5metric *m, double v_r, sim5tetrad *t)
//! Tetrad of observer that moves purely in radial direction.
//! Returns basis vectors of a tetrad e_{(a)}^{\mu} of an observer that moves purely
//! in radial direction (it does have zero phi velocity component). Positive value of `v_r` 
//! means outward motion, a negative value inward motion.
//!
//! Tetrad vector orientation: 
//! e[1] (x-vector) is oriented along increasing r (even when off equatorial plane),
//! e[2] (z-vector) is oriented along decreasing theta (goes "upwards" from eq plane), 
//! e[3] (y-vector) is oriented along increasing phi.
//!
//! @param m metric
//! @param v_r radial velocity component in [c]
//! @param t tetrad
//!
//! @result Returns tetrad vectors in `t`.
{
    if (v_r==0.0) return tetrad_zamo(m, t);

    double g00 = m->g00;
    double g11 = m->g11;
    double U0 = sqrt((-1.-sqr(v_r)*g11)/g00);
    double U1 = v_r;

    t->e[0][0] = U0;
    t->e[0][1] = U1;
    t->e[0][2] = 0.0;
    t->e[0][3] = 0.0;

    double UG = U0*U0*g00 + U1*U1*g11;
    t->e[1][0] = -U1*sqrt(UG*g11*g00)*U0/(g11*UG)*g11/(U0*g00);
    t->e[1][1] = sqrt(UG*g11*g00)*U0/(g11*UG);
    t->e[1][2] = 0.0;
    t->e[1][3] = 0.0;

    t->e[2][0] = 0.0;
    t->e[2][1] = 0.0;
    t->e[2][2] = -1./sqrt(m->g22);  // vector is oriented "upwards" along decreasing theta
    t->e[2][3] = 0.0;

    t->e[3][0] = 0.0;
    t->e[3][1] = 0.0;
    t->e[3][2] = 0.0;
    t->e[3][3] = 1./sqrt(m->g33);

    t->metric = *m;
}



DEVICEFUNC
void tetrad_azimuthal(sim5metric *m, double Omega, sim5tetrad *t)
//! Tetrad of observer that moves purely in azimuthal direction.
//! Returns basis vectors of a tetrad e_{(a)}^{\mu} of an observer that moves purely
//! in azimuthal direction with angular velocity Omega. 
//!
//! Tetrad vector orientation: 
//! e[1] (x-vector) is oriented along increasing r (even when off equatorial plane),
//! e[2] (z-vector) is oriented along decreasing theta (goes "upwards" from eq plane), 
//! e[3] (y-vector) is oriented along increasing phi.
//!
//! @param m metric
//! @param Omega angular velocity in [g.u.]
//! @param t tetrad
//!
//! @result Returns tetrad vectors in `t`.
{
    if (Omega==0.0) return tetrad_zamo(m, t);

    double g00 = m->g00;
    double g33 = m->g33;
    double g03 = m->g03;
    double U0 = sqrt(-1.0/(g00 + 2.*Omega*g03 + sqr(Omega)*g33));
    double U3 = U0*Omega;

    t->e[0][0] = U0;
    t->e[0][1] = 0.0;
    t->e[0][2] = 0.0;
    t->e[0][3] = U3;

    t->e[1][0] = 0.0;
    t->e[1][1] = sqrt(1./m->g11);
    t->e[1][2] = 0.0;
    t->e[1][3] = 0.0;

    t->e[2][0] = 0.0;
    t->e[2][1] = 0.0;
    t->e[2][2] = -sqrt(1./m->g22);  // vector is oriented "upwards" along decreasing theta
    t->e[2][3] = 0.0;

    double k1 = (g03*U3+g00*U0);
    double k2 = (g33*U3+g03*U0);
    t->e[3][0] =  - sign(k1)*k2 / sqrt((g33*g00-g03*g03)*(g00*U0*U0+g33*U3*U3+2.0*g03*U0*U3));
    t->e[3][1] = 0.0;
    t->e[3][2] = 0.0;
    t->e[3][3] = t->e[3][0] * (-k1/k2);

    t->metric = *m;
}



DEVICEFUNC
void tetrad_surface(sim5metric *m, double Omega, double V, double dhdr, sim5tetrad *t)
//! Tetrad of observer that moves along a surface.
//! Returns basis vectors of a tetrad e_{(a)}^{\mu} of an observer that moves along an axisymmetric
//! surface in Kerr spacetime, i.e. the obeserver moves azimuthally with angular velocity Omega and 
//! drifs radially along the surface with velocity V, which itself is measured locally in 
//! the corotating frame. The orientation of the surface is given locally by the derivative dH/dR, 
//! where H is the height of the surface above equatoriual plane [r*cos(theta)] and 
//! R is the BL radial coordinate in equatorial plane [r*sin(theta)].
//!
//! Based on: Sadowski+2011, Appendix A (http://adsabs.harvard.edu/abs/2011A%26A...532A..41S); 
//! note the +--- metric signature used there.
//!
//! Tetrad vector orientation: 
//! e[0] (t-vector) is oriented along the direction of the observer's 4-velocity
//! e[1] (x-vector) is a spacelike vector in [r,theta] plane tangent to the surface and oriented outwards
//! e[2] (z-vector) is a spacelike vector in [r,theta] plane normal to the surface and oriented upwards
//! e[3] (y-vector) is a remaining spacelike vector orthogonal to all the others
//!
//! @param m metric
//! @param Omega angular velocity in [g.u.]
//! @param radial drift velocity mesured in corotating frame [c]
//! @param t tetrad
//!
//! @result Returns tetrad vectors in `t`.

{
    double g00 = m->g00;
    double g11 = m->g11;
    double g22 = m->g22;
    double g33 = m->g33;
    double g03 = m->g03;

    // get contra-variant metric components g^\mu\nu
    sim5metric M;
    kerr_metric_contravariant(m->a, m->r, m->m, &M);


    // zero-radial-velocity surface tangent vector S0:
    // - get the components of space-like surface tangent vector S0 for an observer that
    //   corotates with the fluid (has no radial component of velocity)
    // - S0 is contained in [r,theta] plane and is oriented in the positive radial direction
    // - obviously S = S0*[0,1,d\theta/dr,0], where condition S.S=1 fixes S0
    double S0r = 1.0/sqrt(g11+g22*sqr(dhdr));
    double S0h = S0r*dhdr;

    //fprintf(stderr, "S0.S0 = %e\n", g11*sqr(S0r)+g22*sqr(S0h));

    // compute quantity v (Sadowski+2011, Eq. A.10)
    double ur = V/sqrt(1.-V*V)/sqrt(g11);
    double v = sign(V) * sqrt((sqr(ur/S0r)*(-g00-2.*Omega*g03-sqr(Omega)*g33))/(1.+sqr(ur/S0r)));

    // timelike vector U corresponding to observer's 4-velocity with both azimuthal and radial component
    // U = A(\eta^\mu + Omega*\xi^\mu + v*S_0^\mu); c.f. Sadowski+2011, Eq. A.5
    t->e[0][0] = 1.0;
    t->e[0][1] = v*S0r;
    t->e[0][2] = v*S0h;
    t->e[0][3] = Omega;
    vector_norm_to(t->e[0], -1.0, m); 
    //fprintf(stderr, "U.U = %e\n", dotprod(t->e[0],t->e[0],m));

    // surface tangent vector S
    // this spacelike vector lives in tangent plane to surface and satisfies S.N=0
    // orientation is set such that the vector points in the positive radial direction,
    // i.e. in the limit dhdr=0 it becomes same as ZAMO e[1] vector
    // c.f. Sadowski+2011, Eq. A.12 (signs differ)
    t->e[1][0] = (v*t->e[0][0]);
    t->e[1][1] = (v*t->e[0][1] + S0r/t->e[0][0]);
    t->e[1][2] = (v*t->e[0][2] + S0h/t->e[0][0]);
    t->e[1][3] = (v*t->e[0][3]);
    vector_norm_to(t->e[1], 1.0, m);
    //fprintf(stderr, "S.S = %e\n", dotprod(t->e[1],t->e[1],m));

    // surface normal vector N
    // this spacelike vector lives in [r,theta] plane and satisfies
    //   N^\mu = grad F = [0,\partial{F}/\partial{r},\partial{F}/\partial{theta},0] =
    //   N0*[0, d\theta/dr, -1, 0], where d\theta/dr=-(\partial{F}/\partial{r})/(\partial{F}/\partial{\theta})
    //   and where F(r,theta)=0 defines the surface
    // orientation is set such that the vector points upward,
    // i.e. in the limit dhdr=0 it becomes same as ZAMO e[2] vector
    // c.f. Sadowski+2011, Eq. A.3 (signs differ)
    t->e[2][0] = 0.0;
    t->e[2][1] = dhdr;
    t->e[2][2] = -1.0;
    t->e[2][3] = 0.0;
    vector_norm_to(t->e[2], 1.0, m);
    
    // the remaining [t,phi] plane space-like vector K satisfying K.U=0
    // c.f. Sadowski+2011, Eq. A.8
    t->e[3][0] = -(g03+g33*Omega)/(g00+g03*Omega);
    t->e[3][1] = 0.0;
    t->e[3][2] = 0.0;
    t->e[3][3] = 1.0;
    vector_norm_to(t->e[3], 1.0, m);
    //fprintf(stderr, "K.K = %e\n", dotprod(t->e[3],t->e[3],m));

    //fprintf(stderr, "U.K = %e\n", dotprod(t->e[0],t->e[3],m));
    //fprintf(stderr, "U.N = %e\n", dotprod(t->e[0],t->e[2],m));
    //fprintf(stderr, "U.S = %e\n", dotprod(t->e[0],t->e[1],m));
    //fprintf(stderr, "K.N = %e\n", dotprod(t->e[3],t->e[2],m));
    //fprintf(stderr, "K.S = %e\n", dotprod(t->e[3],t->e[1],m));
    //fprintf(stderr, "N.S = %e\n", dotprod(t->e[2],t->e[1],m));

    t->metric = *m;
}



DEVICEFUNC
void bl2on(double Vin[4], double Vout[4], sim5tetrad* t)
//! Vector transformation from coordinate to local frame.
//! Transforms a vector `V` from coordinate (Boyer-Lindquist) frame to local (orthonormal) 
//! frame specified by tetrad `t`.
//!
//! MATH: V^(a) = e^(a)_\mu * V^\mu, where e^(a)_\mu = e^\nu_(b) * g_\mu\nu * n^ab
//!       V^(a) = dotprod(e_(b)^\mu, Vin^\mu) * n^ab
//!
//! @param Vin vector to transform (in coordinate basis) 
//! @param Vout transformed vector (in local basis)
//! @param t transformation tetrad
//!
//! @result Local vector Vout.
{
    Vout[0] = -dotprod(t->e[0], Vin, &t->metric);
    Vout[1] = +dotprod(t->e[1], Vin, &t->metric);
    Vout[2] = +dotprod(t->e[2], Vin, &t->metric);
    Vout[3] = +dotprod(t->e[3], Vin, &t->metric);
}


DEVICEFUNC
void on2bl(double Vin[4], double Vout[4], sim5tetrad* t)
//! Vector transformation from local to coordinate frame.
//! Transforms a vector `V` from local (orthonormal) frame specified by tetrad `t` to 
//! coordinate (Boyer-Lindquist) frame.
//!
//! MATH: V^\mu = e^\mu_(a) * V^(a)  or V^i = V^j * e_j^i
//! 
//! @param Vin vector to transform (in local basis) 
//! @param Vout transformed vector (in coordinate basis)
//! @param t transformation tetrad
//!
//! @result Coordinate vector Vout.
{
    int i,j;
    for (i=0;i<4;i++) {
        Vout[i] = 0.0;
        for (j=0;j<4;j++) Vout[i] += Vin[j] * t->e[j][i];
    }
    //Vout[0] = + (t->e[0][0]*Vin[0] + t->e[1][0]*Vin[1] + t->e[2][0]*Vin[2] + t->e[3][0]*Vin[3]);
    //Vout[1] = + (t->e[0][1]*Vin[0] + t->e[1][1]*Vin[1] + t->e[2][1]*Vin[2] + t->e[3][1]*Vin[3]);
    //Vout[2] = + (t->e[0][2]*Vin[0] + t->e[1][2]*Vin[1] + t->e[2][2]*Vin[2] + t->e[3][2]*Vin[3]);
    //Vout[3] = + (t->e[0][3]*Vin[0] + t->e[1][3]*Vin[1] + t->e[2][3]*Vin[2] + t->e[3][3]*Vin[3]);
}




//-----------------------------------------------------------------
// orbital motion
//-----------------------------------------------------------------


DEVICEFUNC INLINE
double r_bh(double a)
//! Black hole event horizon radius.
//!
//! @param a black hole spin
//!
//! @result Event horizon radius in [rg].
{
    return 1. + sqrt(1.-sqr(a));
}



DEVICEFUNC INLINE
double r_ms(double a)
//! Radius of marginally stable orbit.
//!
//! @param a black hole spin
//!
//! @result Marginally stable orbit radius in [rg].
{
    double z1 = 1. + sqrt3(1.-sqr(a))*(sqrt3(1.+a) + sqrt3(1.-a));
    double z2 = sqrt(3.*sqr(a) + sqr(z1));
    return 3.+z2-/*+*/sqrt((3.-z1)*(3.+z1+2.*z2));
}


DEVICEFUNC INLINE
double r_mb(double a)
//! Radius of marginally bound orbit.
//! Marginally bound orbit is minimal radius of bound (E<0) circular and parabolic orbits.
//! http://adsabs.harvard.edu/abs/1972ApJ...178..347B, eq. 2.19.
//!
//! @param a black hole spin
//!
//! @result Marginally bound orbit radius in [rg].
{
    return (2.-a) + 2.*sqrt(1.-a);
}


DEVICEFUNC INLINE
double r_ph(double a)
//! Radius of photon orbit.
//! Photon orbit is radius of unstable circular photon orbit.
//! http://adsabs.harvard.edu/abs/1972ApJ...178..347B, eq. 2.18.
//!
//! @param a black hole spin
//!
//! @result Marginally bound orbit radius in [rg].
{
    return 2.0*(1.0+cos(2./3.*acos(-a)));
}



DEVICEFUNC INLINE
double OmegaK(double r, double a)
//! Angular frequency of Keplerian orbital motion.
//!
//! @param r radius [rg]
//! @param a black hole spin
//!
//! @result Angular frequency [Hz].
{
    return 1./(a + pow(r,1.5));
}


DEVICEFUNC INLINE
double ellK(double r, double a)
//! Specific angular momentum of Keplerian orbital motion.
//!
//! @param r radius [rg]
//! @param a black hole spin
//!
//! @result Keplerian specific angular momentum [g.u.].
{
    /*
    double r2 = sqr(r);
    double a2 = sqr(a);
    double D = r2 + a2 - 2.*r;
    double gtt = -1. + 2.0/r;
    double gff = (sqr(r2 + a2) - a2*D)/r2;
    double gtf = -2.*a/r;

    double Omega = OmegaK(r,a);
    return -(gtf + Omega*gff) / (gtt + Omega*gtf);
    */
    
    // formula of Komissarov(2008)
    return (sqr(r)-2.*a*sqrt(r)+sqr(a)) / (sqrt(r)*r-2.*sqrt(r)+a);    
}


DEVICEFUNC INLINE
double omega_r(double r, double a)
//! Angular frequency of radial epyciclic motion.
//!
//! @param r radius [rg]
//! @param a black hole spin
//!
//! @result Angular velocity [Hz].
{
    return OmegaK(r,a) * sqrt(1.-6./r+8.*a/sqrt(r*r*r)-3.*a*a/sqr(r));
}


DEVICEFUNC INLINE
double omega_z(double r, double a)
//! Angular frequency of vertical epyciclic motion.
//!
//! @param r radius [rg]
//! @param a black hole spin
//!
//! @result Angular velocity [Hz].
{
    return OmegaK(r,a) * sqrt(1.-4.*a/sqrt(r*r*r)+3.*a*a/sqr(r));
}


DEVICEFUNC INLINE
double Omega_from_ell(double ell, sim5metric *m)
//! Angular frequency corresponding to given angular momentum.
//!
//! @param ell specific angular momentum [g.u.]
//! @param m metric
//!
//! @result Angular frequency [Hz].
{
    return  -(m->g03 + ell*m->g00) / (m->g33 + ell*m->g03);
}


DEVICEFUNC INLINE
double ell_from_Omega(double Omega, sim5metric *m)
//! Specific angular momentum corresponding to given angular frequency.
//!
//! @param Omega angular frequency [Hz]
//! @param m metric
//!
//! @result Specific angular momentum [g.u.].
{
    return -(m->g03 + m->g33*Omega)/(m->g00 + m->g03*Omega);
}


DEVICEFUNC INLINE
double gfactorK(double r, double a, double l)
//! Redshift factor. 
//! Relativistic correction to energy of photon that is emitted by fluid at Keplerian 
//! rotation in equatorial plane; includes Doppler effect and gravitational redshift.
//!
//! @param r radius [rg]
//! @param a black hole spin
//! @param l photon motion constant lambda
//!
//! @result Redshift factor.
{
    double OmegaK = 1./(a + pow(r,1.5));
    return sqrt(1. - 2./r * pow(1.-a*OmegaK,2.) - (r*r+a*a)*sqr(OmegaK)) / (1. - OmegaK*l);
}



//-----------------------------------------------------------------
// photon motion
//-----------------------------------------------------------------


DEVICEFUNC
void photon_momentum(double a, double r, double m, double l, double q, double r_sign, double m_sign, double k[4])
//! Photon 4-momentum vector .
//! Returns photon 4-momentum vector k^\mu such that k*k=0 (null vector).
//! The orientation of the resulting vector `k` is given by the signs of `r_sign` and `m_sign` parameters.
//! See MTW; Rauch & Blandford (1994), Appendix A; Li+05
//!
//! @param a black hole spin
//! @param r radial coordinate [rg]
//! @param m poloidal coordinate [cos(theta)]
//! @param l photon motion constant lambda
//! @param q photon motion constant Q (Carter's constant)
//! @param r_sign sign of k[1] component of resulting momentum vector 
//! @param m_sign sign of k[2] component of resulting momentum vector 
//! @param k resulting momentum vector (output) 
//!
//! @result Photon momentum vector `k`.
{
    double a2 = sqr(a);
    double l2 = sqr(l);
    double r2 = sqr(r);
    double m2 = sqr(m);
    double S = r2 + a2*m2;
    double D = r2 - 2.*r + a2;

    // after Li+05
    double R = sqr(r2+a2-a*l) - D*( sqr(l-a) + q );        // dr/dl; eq.A2
    double M = q - l2*m2/(1.-m2) + a2*m2;                  // dtheta/dl; eq.A3

    if ((M<0.0) && (-M<1e-8)) M = 0.0;
    if ((R<0.0) && (-R<1e-8)) R = 0.0;

    #ifndef CUDA
    if (R<0.0) error("(photon_momentum): R<0 (R=%.10e)", R);
    if (M<0.0) {
        error("(photon_momentum): M<0 (M=%.10e, l=%.4e  q=%.4e)", M, l, q);
        k[0]=k[1]=k[2]=k[3]=NAN;
        return;
    }
    #endif

    k[0] = +1/S * ( -a*(a*(1.-m2)-l) + (r2+a2)/D*(r2+a2-a*l) );
    k[1] = +1/S * sqrt(R);
    k[2] = +1/S * sqrt(M);
    k[3] = +1/S * ( -a + l/(1.-m2) + a/D*(r2+a2-a*l) );

    if (r_sign<0.0) k[1] = -k[1];
    if (m_sign<0.0) k[2] = -k[2];

    /*
    // eveluate derivative of k if it is requested (dk is not NULL)
    if (dk) {
        int i,b,c;
        double G[4][4][4];
        kerr_connection(a, r, m, G);
        for (i=0;i<4;i++) {
            // dk^\mu = -Gamma^\mu_\alpha\beta k^\alpha k^\beta
            // Gamma is symmetric in \alpha,\beta
            dk[i] = 0.0;
            for (b=0;b<4;b++) for (c=b;c<4;c++) dk[i] += -G[i][b][c]*k[b]*k[c];
        }
    }
    */
}


DEVICEFUNC
void photon_motion_constants(double a, double r, double m, double k[4], double* L, double* Q)
//! Constants of motion L,Q for null geodesic.
//!
//! @param a black hole spin
//! @param r radial coordinate [rg]
//! @param m poloidal coordinate [cos(theta)]
//! @param k photon 4-momentum vector
//! @param L photon motion constant lambda (output)
//! @param Q photon motion constant Q^2 (Carter's constant, output)
//!
//! @result Photon motion constants `L` and `Q`.
{
    double a2 = sqr(a);
    double r2 = sqr(r);
    double s2 = 1.-m*m;
    double D  = r2 - 2.*r + a2;

    double l;
    double nf = k[3]/k[0];
    double nh = sqr(k[2])/sqr(k[0]);

    // Mathematica's solution of "l" for nf==k[3]/k[0]
    *L = l =(-a*a2 + sqr(a2)*nf + nf*sqr(r2) + a*(D-r2) + a2*nf*(2.*r2-D*s2))*s2 /
        (D - a*s2*(a-a2*nf + nf*(D-r2)));

    // Mathematica's solution of "q" for nh==sqr(k[2])/sqr(k[0])
    *Q = pow(a*(l-a*s2) + ((a2+r2)*(a2-a*l+r2))/D, 2.0) *
        (nh - (sqr(D*m)*(sqr(l)-a2*s2))/(-s2*pow(sqr(a2)-a*a2*l+sqr(r2)+a*l*(D-r2)+a2*(2.*r2-D*s2),2.0)));

    #ifndef CUDA
    if (isnan(*L)) warning("ERR (photon_motion_constants): L is NaN (%e, k=%e/%e/%e/%e)\n", *L, k[0], k[1], k[2], k[3]);
    if (isnan(*Q)) warning("ERR (photon_motion_constants): Q is NaN (%e, k=%e/%e/%e/%e)\n", *Q, k[0], k[1], k[2], k[3]);
    #endif
}



DEVICEFUNC
double photon_carter_const(double k[4], sim5metric *metric)
//! Carter's constant Q for null geodesic.
//!
//! @param k photon 4-momentum vector
//! @param m metric
//!
//! @result Carter's constants `Q`.
{
    double m2 = sqr(metric->m);
    double kt = k[0]*metric->g00 + k[3]*metric->g03;
    double kh = k[2]*metric->g22;
    double kf = k[3]*metric->g33 + k[0]*metric->g03;
    return sqr(kh) + sqr(kf)*m2/(1.-m2) - sqr(metric->a)*sqr(kt)*m2;
}




//-----------------------------------------------------------------
// velocity
//-----------------------------------------------------------------


DEVICEFUNC INLINE
void fourvelocity_zamo(sim5metric *m, double U[4])
//! 4-velocity of ZAMO (locally non-rotating) observer.
//!
//! @param m metric
//! @param U 4-velocity (output)
//!
//! @result 4-velocity `U`.
{
    U[0] = sqrt(m->g33/(sqr(m->g03) - m->g33*m->g00));
    U[1] = 0.0;
    U[2] = 0.0;
    U[3] = -U[0] * m->g03/m->g33;
}



DEVICEFUNC INLINE
void fourvelocity_azimuthal(double Omega, sim5metric *m, double U[4])
//! 4-velocity of azimuthally rotating observer.
//!
//! @param Omega angular frequency of circular motion
//! @param m metric
//! @param U 4-velocity (output)
//!
//! @result 4-velocity `U`.
{
    U[0] = sqrt(-1.0/(m->g00 + 2.*Omega*m->g03 + sqr(Omega)*m->g33));
    U[1] = 0.0;
    U[2] = 0.0;
    U[3] = U[0]*Omega;
}



DEVICEFUNC INLINE
void fourvelocity_radial(double vr, sim5metric *m, double U[4])
//! 4-velocity of radially moving observer.
//!
//! @param vr radial component of 4-velocity (dr/dtau) 
//! @param m metric
//! @param U 4-velocity (output)
//!
//! @result 4-velocity `U`.
{
    U[0] = sqrt((-1.0 - sqr(vr)*m->g11)/m->g00);
    U[1] = vr;
    U[2] = 0.0;
    U[3] = 0.0;
}



DEVICEFUNC INLINE
double fourvelocity_norm(
    double U1, double U2, double U3, sim5metric *m)
//*********************************************************
{
    double D = sqr(m->g03*U3) - m->g00*m->g11*sqr(U1) - m->g00*m->g22*sqr(U2) - m->g00*m->g33*sqr(U3) - m->g00;
    return (-m->g03*U3-sqrt(D))/m->g00;
}



DEVICEFUNC 
void fourvelocity(
    double U1, double U2, double U3, sim5metric *m, double U[])
//*********************************************************
{
    double D = sqr(m->g03*U3) - m->g00*m->g11*sqr(U1) - m->g00*m->g22*sqr(U2) - m->g00*m->g33*sqr(U3) - m->g00;
    double N = (-m->g03*U3 - sqrt(D))/m->g00;
    U[0] = 1./N;
    U[1] = U1/N;
    U[2] = U2/N;
    U[3] = U3/N;
}








/* 
replaced by tetrad_general

DEVICEFUNC
void ortho_tetrad_U(
    double U[4],
    double g00, double g11, double g22, double g33, double g03,
    double e0[4], double e1[4], double e2[4], double e3[4])
// calculates a convenient orthonormal tetrad for a fluid with
// four-velocity U, which satisfies condition U^\phi >> U^r, U^\theta
// see Krolik, Hawley & Hirose (2005), ApJ, 622, 1088
//     Beckwith, Hawley, Krolik (2008), MNRAS, 390, 21 (corrected formulae)
// - puts y-vector oriented along increasing phi
{

    double k1 = sqrt(fabs(g33*sqr(U[3]) + g11*sqr(U[1]) + U[0]*(2.*g03*U[3]+g00*U[0])));
    double k2 = sqrt(fabs(g33*sqr(U[3]) + g11*sqr(U[1]) + g22*sqr(U[2]) + U[0]*(2.*g03*U[3]+g00*U[0])));
    double k3 = sqrt(fabs(g33*sqr(U[3]) + U[0]*(2.*g03*U[3]+g00*U[0])));
    double f;

    e0[0] = U[0];
    e0[1] = U[1];
    e0[2] = U[2];
    e0[3] = U[3];

    f = +1./(k1*k3);
    e1[0] = f * sqrt(g11)*U[1]*U[0];
    e1[1] = f * sqr(k3)/sqrt(g11);
    e1[2] = 0.0;
    e1[3] = f * sqrt(g11)*U[1]*U[3];

    f = -1./(k1*k2);
    e2[0] = f * sqrt(g22)*U[2]*U[0];
    e2[1] = f * sqrt(g22)*U[2]*U[1];
    e2[2] = f * sqr(k1)/sqrt(g22);
    e2[3] = f * sqrt(g22)*U[2]*U[3];

    f = +1./(k3*sqrt(fabs(g33*g00-g03*g03)));
    e3[0] = f * (+g33*U[3] + g03*U[0]);
    e3[1] = 0.0;
    e3[2] = 0.0;
    e3[3] = f * (-g03*U[3] - g00*U[0]);
}
*/


/*
commented out because orientation of vecotrs is unclear

DEVICEFUNC
void ortho_tetrad_U_phi_r_motion(
    double U[4],
    double g00, double g11, double g22, double g33, double g03,
    double e0[4], double e1[4], double e2[4], double e3[4])

// calculates a convenient orthonormal tetrad for a fluid with
// four-velocity U, which satisfies condition U^\theta = 0
// - puts x-vector oriented along increasing r, but tilted into phi by U^r boost
// - puts z-vector oriented along decreasing theta (goes "upwards" from eq plane)
// - puts y-vector oriented along increasing phi
{

    double k1 = sqrt(fabs(g33*sqr(U[3]) + g11*sqr(U[1]) + U[0]*(2.*g03*U[3]+g00*U[0])));
    double k3 = sqrt(fabs(g33*sqr(U[3])                 + U[0]*(2.*g03*U[3]+g00*U[0])));
    double f;

    e0[0] = U[0];
    e0[1] = U[1];
    e0[2] = 0.0;
    e0[3] = U[3];

    f = +1./(k1*k3);
    e1[0] = f * sqrt(g11)*U[1]*U[0];
    e1[1] = f * sqr(k3)/sqrt(g11);
    e1[2] = 0.0;
    e1[3] = f * sqrt(g11)*U[1]*U[3];

    e2[0] = 0.0;
    e2[1] = 0.0;
    e2[2] = -1./sqrt(g22);
    e2[3] = 0.0;

    f = +1./(k3*sqrt(fabs(g33*g00-g03*g03)));
    e3[0] = f * (+g33*U[3] + g03*U[0]);
    e3[1] = 0.0;
    e3[2] = 0.0;
    e3[3] = f * (-g03*U[3] - g00*U[0]);
}
*/



