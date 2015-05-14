//************************************************************************
//    SIM5 library
//    sim5kerr.c - basic properties of Kerr spacetime
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************

/*
#include "sim5config.h"
#ifndef CUDA
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim5math.h"
#include "sim5kerr.h"
#endif
*/


//-----------------------------------------------------------------
// metric, tetrads and vectors
//-----------------------------------------------------------------

static long prof_N = 0;
static double prof_t = 0.0;

void sim5kerr_profile() {
    fprintf(stderr, "sim5kerr_profile: N=%ld t=%.2e t/N=%.3e\n", prof_N, prof_t, prof_t/(double)(prof_N));
}



DEVICEFUNC
void flat_metric(double a, double r, double m, sim5metric *metric)
//*********************************************************
// returns covariant Kerr metric g_\mu\nu
// inputs: spin <a>, radius <r>, cos(theta) <m>
// output: Kerr metric <metric>
{
    metric->a = a;
    metric->r = r;
    metric->m = m;
    metric->g00 = -1;
    metric->g11 = 1;
    metric->g22 = r*r;
    metric->g33 = r*r*(1.-m*m); 
    metric->g03 = 0.0;
}



DEVICEFUNC
void kerr_metric(double a, double r, double m, sim5metric *metric)
//*********************************************************
// returns covariant Kerr metric g_\mu\nu
// inputs: spin <a>, radius <r>, cos(theta) <m>
// output: Kerr metric <metric>
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
void kerr_connection(double a, double r, double m, double G[4][4][4])
//*********************************************************
// returns Christoffel symbols for Kerr metric connection (Gamma^\mu_\alpha\beta)
// inputs: spin <a>, radius <r>, cos(theta) <m>
// output: connection coefficient matrix
// NOTE: Christoffel tensor Gamma^i_(jk) is symmetric in lower two indices (j,k).
//       (!) This function only evaluates components of the tensor, where j<k and 
//       multiplies the value of these coeficients by a factor of two.
//       In this way, both summation over j=1..4,k=1..4 and over j=1..4,k=1..j
//       will give the same result.
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
void flat_connection(double a, double r, double m, double G[4][4][4])
//*********************************************************
// returns Christoffel symbols for metric connection (Gamma^\mu_\alpha\beta)
// in the limit of M=0 and a=0 (in flat spacetime with spherical coordinates)
// inputs: spin <a>, radius <r>, cos(theta) <m>
// output: connection coefficient matrix
// NOTE: (!) the same note as for kerr_connection() function applies
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
   

DEVICEFUNC INLINE
void Gamma(double G[4][4][4], double v[4], double result[4]) {
    int i,j,k;
    for (i=0;i<4;i++) {
        result[i] = 0.0;
        for (j=0;j<4;j++) for (k=j;k<4;k++) result[i] += -G[i][j][k]*v[j]*v[k];
    }
}


DEVICEFUNC INLINE
void vector_covariant(double V1[4], double V2[4], sim5metric* m) 
//*********************************************************
// returns covariant version of a vector
// inputs: contravariant vector <V1> (k^\mu), metric <m>
// output: covariant vector <V2> (k_\mu)
// NOTE: if metric is NULL, flat space metric is used
{
    V2[0] = V1[0]*m->g00 + V1[3]*m->g03;
    V2[1] = V1[1]*m->g11;
    V2[2] = V1[2]*m->g22;
    V2[3] = V1[3]*m->g33 + V1[0]*m->g03;
}


DEVICEFUNC INLINE
double vector_norm(double V[4], sim5metric* m) 
//*********************************************************
// returns vector norm: sqrt(V*V)
// inputs: vector <V1> (k^\mu), metric <m>
// output: vector norm
// NOTE: if metric is NULL, flat space metric is used
// NOTE: only works for space-like vectors, where V*V>0
{
    return sqrt(dotprod(V, V, m));
}



DEVICEFUNC INLINE
void vector_norm_to(double V[4], sim5metric* m, double norm)
//*********************************************************
{
    double N = dotprod(V, V, m);
    V[0] *= sqrt(fabs(norm/N));
    V[1] *= sqrt(fabs(norm/N));
    V[2] *= sqrt(fabs(norm/N));
    V[3] *= sqrt(fabs(norm/N));
}


DEVICEFUNC INLINE
double dotprod(double V1[4], double V2[4], sim5metric* m) 
//*********************************************************
// returns scalar product of two vectors
// inputs: vector <V1>, vector <V2>, metric <m>
// output: result of the scalar product
// NOTE: if metric is NULL, flat space metric is used
{
    if (m) 
        return V1[0]*V2[0]*m->g00 + V1[1]*V2[1]*m->g11 + V1[2]*V2[2]*m->g22 + 
               V1[3]*V2[3]*m->g33 + V1[0]*V2[3]*m->g03 + V1[3]*V2[0]*m->g03;
    else
        return -V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2] + V1[3]*V2[3];
}



DEVICEFUNC
void tetrad_zamo(sim5metric *m, sim5tetrad *t)
//*********************************************************
// returns basis vectors for ZAMO observer tetrad e_{(a)}^{\mu}
// inputs: spin <a>, radius <r>, cos(theta) <m>
// output: tetrad <t>
// NOTE: x-vector is oriented along increasing r (even when off equatorial plane)
//       z-vector is oriented along decreasing theta (goes "upwards" from eq plane)
//       y-vector oriented along increasing phi
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
    
    t->m = *m;
}



DEVICEFUNC
void tetrad_radial(sim5metric *m, double v_r, sim5tetrad *t)
//*********************************************************
// returns basis vectors for tetrad e_{(a)}^{\mu} of an observer moving in 
// radial direction with some velocity
// inputs: metric <m>, velocity <v_r>
// output: tetrad <t>
// NOTE: x-vector is oriented along increasing r (even when off equatorial plane)
//       z-vector is oriented along decreasing theta (goes "upwards" from eq plane)
//       y-vector oriented along increasing phi
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
    
    t->m = *m;
}



DEVICEFUNC
void tetrad_azimuthal(sim5metric *m, double Omega, sim5tetrad *t)
//*********************************************************
// returns basis vectors for tetrad e_{(a)}^{\mu} of an observer rotating 
// azimuthally with angular velocity Omaga
// inputs: metric <m>, ang.velocity <Omega>
// output: tetrad <t>
// NOTE: x-vector is oriented along increasing r (even when off equatorial plane)
//       z-vector is oriented along decreasing theta (goes "upwards" from eq plane)
//       y-vector oriented along increasing phi
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
    
    t->m = *m;
}



DEVICEFUNC
void tetrad_surface(sim5metric *m, double Omega, double V, double dhdr, sim5tetrad *t)
//*********************************************************
// returns basis vectors for tetrad e_{(a)}^{\mu} of an observer moving along an axi-symmetric surface in Kerr spacetime
// based on: Sadowski+2011, Appendix A (http://adsabs.harvard.edu/abs/2011A%26A...532A..41S); note the +--- metric signature used there
// the obeserver moves azimuthally with angular velocity Omega and drifs radially with velocity V (measured in CRF frame)
// inputs: metric <m>, ang.velocity <Omega>, radial velocity <V>, surface derivative <dhdr>=d(theta)/d(r)
// output: tetrad <t>
// NOTE: t-vector is oriented along the direction of the observer's 4-velocity
//       x-vector is a spacelike vector in [r,theta] plane tangent to the surface and oriented outwards
//       y-vector is a spacelike vector in [r,theta] plane normal to the surface and oriented upwards
//       z-vector a remaining spacelike vector orthogonal to all the others
{
    double g00 = m->g00;
    double g11 = m->g11;
    double g22 = m->g22;
    double g33 = m->g33;
    double g03 = m->g03;

    // zero-radial-velocity surface tangent vector S0:
    // - get the components of space-like surface tangent vector S0 for an observer that 
    //   corotates with the fluid (has no radial component of velocity)
    // - S0 is contained in [r,theta] plane and satisfies condition S.N=0
    // - orientation is set such that the vector points in the positive radial direction
    // c.f. Sadowski+2011, Eq. A.4 (signs differ)
    double S0r = 1./sqrt(g11 + sqr(g11)/g22*sqr(dhdr));
    double S0h = g11/g22*dhdr*S0r;
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
    vector_norm_to(t->e[0], m, 1.0);  // note: the actual computed norm will be -1.0 because of time-likeness
    //fprintf(stderr, "U.U = %e\n", dotprod(t->e[0],t->e[0],m));

    // surface tangent vector S
    // this spacelike vector lives in [r,theta] plane and satisfies S.N=0
    // orientation is set such that the vector points in the positive radial direction,
    // i.e. in the limit dhdr=0 it becomes same as ZAMO e[1] vector
    // c.f. Sadowski+2011, Eq. A.12 (signs differ)
    t->e[1][0] = (v*t->e[0][0]); 
    t->e[1][1] = (v*t->e[0][1] + S0r/t->e[0][0]); 
    t->e[1][2] = (v*t->e[0][2] + S0h/t->e[0][0]); 
    t->e[1][3] = (v*t->e[0][3]); 
    vector_norm_to(t->e[1], m, 1.0);
    //fprintf(stderr, "S.S = %e\n", dotprod(t->e[1],t->e[1],m));
    
    // surface normal vector N
    // this spacelike vector lives in [r,theta] plane as satisfies N=[0,dF/dr,dF/dtheta,0], where F(r,theta)=0 defines the surface
    // orientation is set such that the vector points upward,
    // i.e. in the limit dhdr=0 it becomes same as ZAMO e[2] vector
    // c.f. Sadowski+2011, Eq. A.3 (signs differ)
    t->e[2][0] = 0.0; 
    t->e[2][1] = dhdr;
    t->e[2][2] = -1.0; 
    t->e[2][3] = 0.0; 
    vector_norm_to(t->e[2], m, 1.0);
    //fprintf(stderr, "N.N = %e\n", dotprod(t->e[2],t->e[2],m));
    
    // the remaining [t,phi] plane space-like vector K satisfying K.U=0
    // c.f. Sadowski+2011, Eq. A.8
    t->e[3][0] = -(g03+g33*Omega)/(g00+g03*Omega);
    t->e[3][1] = 0.0; 
    t->e[3][2] = 0.0; 
    t->e[3][3] = 1.0;
    vector_norm_to(t->e[3], m, 1.0);
    //fprintf(stderr, "K.K = %e\n", dotprod(t->e[3],t->e[3],m));

    //fprintf(stderr, "U.K = %e\n", dotprod(t->e[0],t->e[3],m));
    //fprintf(stderr, "U.N = %e\n", dotprod(t->e[0],t->e[2],m));
    //fprintf(stderr, "U.S = %e\n", dotprod(t->e[0],t->e[1],m));
    //fprintf(stderr, "K.N = %e\n", dotprod(t->e[3],t->e[2],m));
    //fprintf(stderr, "K.S = %e\n", dotprod(t->e[3],t->e[1],m));
    //fprintf(stderr, "N.S = %e\n", dotprod(t->e[2],t->e[1],m));
    
    t->m = *m;
}



DEVICEFUNC
void bl2on(double Vin[4], double Vout[4], sim5tetrad* t)
//*********************************************************
// transforms 4-vector from coordinate to local frame
// inputs: 4-vector components in coordinate frame tetrad <Vin>
// output: 4-vector components in local frame tetrad <Vout>
// MATH: V^(a) = e^(a)_\mu * V^\mu, where e^(a)_\mu = e^\nu_(b) * g_\mu\nu * n^ab
//       V^(a) = dotprod(e_(b)^\mu, Vin^\mu) * n^ab
{
    Vout[0] = -dotprod(t->e[0], Vin, &t->m);
    Vout[1] = +dotprod(t->e[1], Vin, &t->m);
    Vout[2] = +dotprod(t->e[2], Vin, &t->m);
    Vout[3] = +dotprod(t->e[3], Vin, &t->m);
}


DEVICEFUNC
void on2bl(double Vin[4], double Vout[4], sim5tetrad* t)
//*********************************************************
// transforms 4-vector from local frame to coordinate frame
// inputs: 4-vector components in local frame tetrad <Vin>
// output: 4-vector components in coordinate frame tetrad <Vout>
// MATH: V^\mu = e^\mu_(a) * V^(a)  or V^i = V^j * e_j^i
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
//*********************************************************
// event horizon
// inputs: spin <a>
// output: radius of r_ms
{
    return 1. + sqrt(1.-sqr(a));
}



DEVICEFUNC INLINE 
double r_ms(double a)
//*********************************************************
// marginally stable orbit (isco)
// inputs: spin <a>
// output: radius of r_ms
{
    double z1 = 1. + sqrt3(1.-sqr(a))*(sqrt3(1.+a) + sqrt3(1.-a));
    double z2 = sqrt(3.*sqr(a) + sqr(z1));
    return 3.+z2-/*+*/sqrt((3.-z1)*(3.+z1+2.*z2));
}


DEVICEFUNC INLINE 
double r_mb(double a)
//*********************************************************
// marginally bound orbit (minimal radius of bound (E<0) circular and parabolic orbits)
// inputs: spin <a>
// output: radius of r_mb
// http://adsabs.harvard.edu/abs/1972ApJ...178..347B, eq. 2.19
{
    return (2.-a) + 2.*sqrt(1.-a);
}


DEVICEFUNC INLINE 
double r_ph(double a)
//*********************************************************
// photon orbit
// inputs: spin <a>
// output: radius of r_ph
// http://adsabs.harvard.edu/abs/1972ApJ...178..347B, eq. 2.18
{
    return 2.0*(1.0+cos(2./3.*acos(-a)));
}



DEVICEFUNC INLINE 
double OmegaK(double r, double a)
{
    return 1./(a + pow(r,1.5));
}


DEVICEFUNC INLINE
double ellK(double r, double a)
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
    return (sqr(r)-2.*a*sqrt(r)+sqr(a)) / (sqrt(r)*r-2.*sqrt(r)+a);    // Komissarov(2008)
}


DEVICEFUNC INLINE 
double omega_r(double r, double a)
{
    return OmegaK(r,a) * sqrt(1.-6./r+8.*a/sqrt(r*r*r)-3.*a*a/sqr(r));
}


DEVICEFUNC INLINE 
double omega_z(double r, double a)
{
    return OmegaK(r,a) * sqrt(1.-4.*a/sqrt(r*r*r)+3.*a*a/sqr(r));
}


DEVICEFUNC INLINE 
double Omega_from_ell(double ell, sim5metric *m)
{
    return  -(m->g03 + ell*m->g00) / (m->g33 + ell*m->g03);
}


DEVICEFUNC INLINE 
double ell_from_Omega(double Omega, sim5metric *m)
{
    return -(m->g03 + m->g33*Omega)/(m->g00 + m->g03*Omega);
}





//-----------------------------------------------------------------
// photon motion
//-----------------------------------------------------------------


DEVICEFUNC
void photon_momentum(double a, double r, double m, double l, double q2, double r_sign, double m_sign, double k[4])
//*********************************************************
// returns photon 4-momentum vector k^\mu
// inputs: spin <a>, radius <r>, cos(theta) <m>, motion constats <l,q2>, signs
// output: 4-momentum vector <k> (k^\mu such that k*k=0)
//         derivative of k <dk> (dk^\mu / d\lambda)
//         <dk> is optional and can be set to NULL on input in which case it is not evaluated and returned
// see MTW; Rauch & Blandford (1994), Appendix A; Li+05
// NOTE: the orientation of the vector is given by the sign of <r_sign> and <m_sign> parameters
{
    double a2 = sqr(a);
    double l2 = sqr(l);
    double r2 = sqr(r);
    double m2 = sqr(m);
    double S = r2 + a2*m2;
    double D = r2 - 2.*r + a2;

    // after Li+05
    double R = sqr(r2+a2-a*l) - D*( sqr(l-a) + q2 );       // dr/dl; eq.A2 (with q2 = Q)
    double M = q2 - l2*m2/(1.-m2) + a2*m2;                 // dtheta/dl; eq.A3

    if ((M<0.0) && (-M<1e-8)) M = 0.0;
    if ((R<0.0) && (-R<1e-8)) R = 0.0;
    if (R<0.0) fprintf(stderr, "ERR (photon_momentum): R<0 (%.10e)\n", R);
    if (M<0.0) fprintf(stderr, "ERR (photon_momentum): M<0 (%.10e)\n", M);

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
//*********************************************************
// returns constants of motion L,Q for null geodesic
// inputs: spin <a>, radius <r>, cos(theta) <m>, 4-momentum vector <k>
// output: motion constants lambda (L) and q^2 (Q, Carter constant)
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
    
    // Mathematica's solution of "q2" for nh==sqr(k[2])/sqr(k[0])
    *Q = pow(a*(l-a*s2) + ((a2+r2)*(a2-a*l+r2))/D, 2.0) *
        (nh - (sqr(D*m)*(sqr(l)-a2*s2))/(-s2*pow(sqr(a2)-a*a2*l+sqr(r2)+a*l*(D-r2)+a2*(2.*r2-D*s2),2.0)));
    if (isnan(*L)) fprintf(stderr,"#L=%e (%e/%e/%e)\n",*L,k[1],k[2],k[3]);
    if (isnan(*Q)) fprintf(stderr,"#Q=%e (%e/%e/%e)\n",*Q,k[1],k[2],k[3]);
}



DEVICEFUNC
double photon_carter(sim5metric *metric, double k[4])
//*********************************************************
// returns Carter constants for a null geodesic 
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
void fourvelocity_azimuthal(double Omega, sim5metric *m, double U[4])
//*********************************************************
{
    U[0] = sqrt(-1.0/(m->g00 + 2.*Omega*m->g03 + sqr(Omega)*m->g33));
    U[1] = 0.0;
    U[2] = 0.0;
    U[3] = U[0]*Omega;
}



DEVICEFUNC INLINE
void fourvelocity_radial(double vr, sim5metric *m, double U[4])
//*********************************************************
{
    U[0] = sqrt((-1.0 - sqr(vr)*m->g11)/m->g00);
    U[1] = vr;
    U[2] = 0.0;
    U[3] = 0.0;
}























DEVICEFUNC INLINE 
void bl(double p[4], double t, double r, double m, double phi)
{
    p[0] = t;
    p[1] = r;
    p[2] = m;
    p[3] = phi;
}




DEVICEFUNC
void kerr_metric2(double a, double r, double m, double *gTT, double *gRR, double *gHH, double *gFF, double *gTF)
//*********************************************************
{
    double r2 = sqr(r);
    double a2 = sqr(a);
    double m2 = sqr(m);
    double S = r2 + a2*m2;
    double D = r2 - 2.*r + a2;
    double A = sqr(r2 + a2) - a2*D*(1.-m2);
    *gTT = -A/(S*D);
    *gRR = D/S;
    *gHH = 1/S;
    *gFF = (D - a2*(1.-m2))/(D*S*(1.-m2)); 
    *gTF = -2.*a*r/(D*S);
}


DEVICEFUNC
void ortho_tetrad_U(
    double U[4],
    double g00, double g11, double g22, double g33, double g03,
    double e0[4], double e1[4], double e2[4], double e3[4])
//*********************************************************
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

DEVICEFUNC
void ortho_tetrad_U_phi_r_motion(
    double U[4],
    double g00, double g11, double g22, double g33, double g03,
    double e0[4], double e1[4], double e2[4], double e3[4])
//*********************************************************
// calculates a convenient orthonormal tetrad for a fluid with 
// four-velocity U, which satisfies condition U^\theta = 0
// - puts x-vector oriented along increasing r, but tilted into phi by U^r boost
// - puts z-vector oriented along decreasing theta (goes "upwards" from eq plane)
// - puts y-vector oriented along increasing phi
{
// orientation of vecotrs is unclear
/*
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
*/
}


DEVICEFUNC INLINE
double fourvelocity_norm(
    double U1, double U2, double U3,
    double g00, double g11, double g22, double g33, double g03)
//*********************************************************
{
    double D = sqr(g03*U3) - g00*g11*sqr(U1) - g00*g22*sqr(U2) - g00*g33*sqr(U3) - g00;
    return (-g03*U3-sqrt(D))/g00;
}


// ------- polarization routines --------------



DEVICEFUNC
void kappa_pw(double a, double r, double m, double k[4], double f[4], double *kappa1, double *kappa2)
// following Connors, Piran, Stark (1980): kappa_pw = kappa1 - I*kappa2 = (A1 - I*A2)*(r - I*a*cos(theta))
// => kappa1 = +r*A1 - a*cos(theta)*A2; kappa2 = -r*A2 - a*cos(theta)*A1
// returns kappa_pw = kappa1 + I*kappa2
// !! note the definition of kappa1 & kappa2, which os opposite to CPS(1980)
{
    // following Connors, Piran, Stark (1980):
    double A1 = (k[0]*f[1]-k[1]*f[0]) + a*(1.-m*m)*(k[1]*f[3]-k[3]*f[1]);
    double A2 = sqrt(1.-m*m)*( (r*r+a*a)*(k[3]*f[2]-k[2]*f[3]) - a*(k[0]*f[2]-k[2]*f[0]) );
    *kappa1 = +r*A1 - a*m*A2;
    *kappa2 = -r*A2 - a*m*A1;
}


DEVICEFUNC
void stokes_infty(double a, double inc, double alpha, double beta, double kappa1, double kappa2, double *pol_angle)
{
// following Connors, Piran, Stark (1980): (note the opposite definition of kappa1,kappa2)
//    double S = l/sin(inc) - a*sin(inc) = -alpha - a*sin(inc)
//    double T = sqrt(q - pow(l*cos(inc)/sin(inc),2.) + pow(a*cos(inc),2.)) = beta
//    this should be considered for k_\theta<0: if (k_\theta|_\infty < 0) T = -T
    double S = -alpha - a*sin(inc); 
    double T = +beta;
    double X = (-S*kappa2 - T*kappa1)/(S*S+T*T);
    double Y = (-S*kappa1 + T*kappa2)/(S*S+T*T);
    (*pol_angle) = atan2(Y,X);
}


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
