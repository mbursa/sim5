/*
#include "sim5config.h"
#ifndef CUDA
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim5math.h"
#include "sim5utils.h"
#include "sim5polyroots.h"
#include "sim5elliptic.h"
#include "sim5kerr.h"
#include "sim5kerr-geod.h"
#endif
*/


#ifdef CUDA

__host__ __device__ void error(char *s) { 
    //abort(); 
}

#endif


//---------------------------------------------------------------------------
// geodesic rutines -- source-to-infinity 
//---------------------------------------------------------------------------


DEVICEFUNC
void geodesic_s2i_setup(double i, double a, double alpha, double beta, geodesic *g, int *status)
{
    *status = 0;
    if (a < 1e-8) a = 1e-8;

    g->a = a;
    g->i = i;
    g->cos_i = cos(i);
    g->alpha = alpha;
    g->beta  = beta;
    
    g->l = -alpha*sin(i);
    g->q = sqr(beta) + sqr(cos(i))*(sqr(alpha)-sqr(a));
    
    if (g->q == 0.0) {
        // we are not yet ready to handle this case (see Agol&Dexter for solution)
        g->nrr = -1;
        return;
    }

    double q  = g->q;
    double l  = g->l;
    double a2 = sqr(g->a);
    double l2 = sqr(g->l);

    // R integral
    double C = sqr(a-l)+q;
    double D = 2./3.*(q+l2-a2);
    double E = 9./4.*sqr(D)-12.*a2*q;
    double F = -27./4.*sqr3(D)-108.*a2*q*D+108.*sqr(C);
    double X = sqr(F)-4.*sqr3(E);
    double A;
    if (X >= 0) {
        A = (F>sqrt(X)?+1:-1)*1./3.*pow(fabs(F-sqrt(X))/2.,1./3.)+ (F>-sqrt(X)?+1:-1)*1./3.*pow(fabs(F+sqrt(X))/2.,1./3.);
    } else {
         double Z = sqrt(pow(F/54.,2) + pow(sqrt(-X)/54.,2));
         double z = atan2(sqrt(-X)/54.,F/54.);
         A = pow(Z,1./3.)*2.*cos(z/3.);
    }
    double B = sqrt(A+D);
    g->r1 = +B/2. + .5*csqrt(makeComplex(-A+2.*D-4.*C/B,0.0)); 
    g->r2 = +B/2. - .5*csqrt(makeComplex(-A+2.*D-4.*C/B,0.0));
    g->r3 = -B/2. + .5*csqrt(makeComplex(-A+2.*D+4.*C/B,0.0));
    g->r4 = -B/2. - .5*csqrt(makeComplex(-A+2.*D+4.*C/B,0.0));

    sort_roots(&g->nrr, &g->r1, &g->r2, &g->r3, &g->r4);

    if (g->nrr == 0) {
        // we do not yet know how to handle the case of non-equatorial photons
        // skip it silently
        #ifndef CUDA
        //fprintf(stderr, "ERR: all R-integral roots are complex (geodesic_s2i_setup)\n");
        #endif
        return;
    }
    
    
    // find roots of the T-intergal - M = q + (a2 - l2 - q)m2 - a2*m4 = -a2*(m4 + m2*(q + l2 -a2)/a2 - q/a2) = a2(m2m+m2)(m2p-m2)
    // the following straightforward calculation fails for small spins due to numerical cancelation,
    // so that it must be calculated with extra care and the help of equality: m2m*m2p=q/a^2 
    //   g->m2m = 1./(2.*a2) * ( sqrt(sqr(qla)+4.*q*a2) + qla );
    //   g->m2p = 1./(2.*a2) * ( sqrt(sqr(qla)+4.*q*a2) - qla );
    // note: Agol&Dexter define m2m with opposite sign
    double qla = q + l2 - a2;
    X = sqrt(sqr(qla)+4.*q*a2) + qla;
    g->m2m = X/2./a2;
    g->m2p = 2.*q/X;

    if (q > 0.0) {   // photons crossing eq plane
        // both m2p and m2m must be non-negative (since m2p*m2m=q/a^2>0)
        if (!ensure_range(&g->m2p, 0.0, 1.0, 1e-5)) {
            #ifndef CUDA
            fprintf(stderr,"ERR: Q>0 && m2p out of range [0..1] (m2p=%e m2m=%e l=%e q=%e alpha=%e beta=%e)\n", g->m2p, g->m2m, l, q, alpha, beta); 
            #endif
            return;
        }
    }

    if (q < 0.0) {   // photons not crossing eq plane
    }

    // for q>0 (photons crossing eq plane), both m2p and m2m are non-negative (m2p*m2m=q/a^2)
    if (!ensure_range(&g->m2p, 0.0, 1.0, 1e-5)) {
        #ifndef CUDA
        fprintf(stderr,"ERR:m2p<0 (%e/%f/%f)\n", g->m2p,alpha,beta); 
        #endif
        return;
    }

    if (g->m2m >= 0.0) {
        g->m2 = g->m2p/(g->m2p+g->m2m);
//        fprintf(stderr, "m2=%e (q=%e)\n", g->m2,q);
        if (g->m2>=1.0) {*status=0; fprintf(stderr,"WRN:g->m2>=1.0\n");return;}
        if (!ensure_range(&g->m2, 0.0, 1.0, 1e-5)) {
            #ifndef CUDA
            //fprintf(stderr,"m2<0 (%e)\n", g->m2); 
            #endif
            return;
        }
    } else {
        g->m2 = g->m2p/(g->m2p-g->m2m);
        //g->a = 1e-3; a2 = 1e-6;
        //g->q = q = sqr(beta) + sqr(cos(i))*(sqr(alpha)-sqr(g->a));
        #ifndef CUDA
        //fprintf(stderr,"m2m<0 (%e)\n", g->m2m); 
        //fprintf(stderr,"m2m=%.4e   m2p=%.4e   m2=%.4e   qla=%.4e   X=%.4e\n", g->m2m, g->m2p, g->m2, qla, X); 
        #endif
        return;
    }

    g->mK = 1./sqrt(a2*(g->m2p+g->m2m));

    // periastron
    g->rp   = (creal(g->r1) > 0.0) ? creal(g->r1) : 0.0;
    g->rp_R = (creal(g->r1) > 0.0) ? geodesic_s2i_int_R_periastron(g) : 0.0;

    *status = 1;
}


DEVICEFUNC
void geodesic_s2i_setup_local(double a, double l, double q, geodesic *g, int *status)
{
    *status = 1;

    g->a = a;
    g->i = 0.0;
    g->cos_i = 0.0;
    g->alpha = 0.0;
    g->beta  = 0.0;
    
    g->l = l;
    g->q = q;
    
    double a2 = sqr(g->a);
    double l2 = sqr(g->l);


    // find roots of the R-integral
    double C = sqr(a-l)+q;
    double D = 2./3.*(q+l2-a2);
    double E = 9./4.*sqr(D)-12.*a2*q;
    double F = -27./4.*sqr3(D)-108.*a2*q*D+108.*sqr(C);
    double X = sqr(F)-4.*sqr3(E);
    double A;
    if (X >= 0) {
        A = (F>sqrt(X)?+1:-1)*1./3.*pow(fabs(F-sqrt(X))/2.,1./3.)+ (F>-sqrt(X)?+1:-1)*1./3.*pow(fabs(F+sqrt(X))/2.,1./3.);
    } else {
         double Z = sqrt(pow(F/54.,2) + pow(sqrt(-X)/54.,2));
         double z = atan2(sqrt(-X)/54.,F/54.);
         A = pow(Z,1./3.)*2.*cos(z/3.);
    }
    double B = sqrt(A+D);
    g->r1 = +B/2. + .5*csqrt(makeComplex(-A+2.*D-4.*C/B,0.0)); 
    g->r2 = +B/2. - .5*csqrt(makeComplex(-A+2.*D-4.*C/B,0.0));
    g->r3 = -B/2. + .5*csqrt(makeComplex(-A+2.*D+4.*C/B,0.0));
    g->r4 = -B/2. - .5*csqrt(makeComplex(-A+2.*D+4.*C/B,0.0));
    sort_roots(&g->nrr, &g->r1, &g->r2, &g->r3, &g->r4);
    
    if (g->nrr == 0) {
        *status = 0;
        #ifndef CUDA
        // we do not yet know how to handle the case of non-equatorial photons
        // skip it silently
        //fprintf(stderr, "ERR: all R-integral roots are complex (geodesic_s2i_setup)\n");
        #endif
        return;
    }


    //find roots of the T-intergal
    // extra care must be taken here for small spins
    double qla = q + l2 - a2;
//???    if (qla<0.0) {fprintf(stderr,"qla2<0 a=%.2e, b=%.2e\n", alpha,beta); *status=0; return;}
    // the following straightforward calculation fails for small spins due to numerical cancelation,
    // so that it must be calculated with extra case and the help of equality: m2m*2mp=q/a^2 
    // g->m2m = 1./(2.*a2) * ( sqrt(sqr(qla)+4.*q*a2) + qla );
    // g->m2p = 1./(2.*a2) * ( sqrt(sqr(qla)+4.*q*a2) - qla );
    X = sqrt(sqr(qla)+4.*q*a2) + qla;
    g->m2m = X/2./a2;
    g->m2p = 2.*q/X;
    // for q>0 (photons crossing eq plane), both m2p and m2m are non-negative (m2p*m2m=q/a^2)
    if (!ensure_range(&g->m2p, 0.0, 1.0, 1e-5)) {
        #ifndef CUDA
        fprintf(stderr,"m2p<0 (%e/%f/%f)\n", g->m2p,l,q); 
        #endif
        *status=0; 
        return;
    }

    g->m2 = g->m2p/(g->m2p+g->m2m);
    if ((g->q>0)&&(!ensure_range(&g->m2, 0.0, 1.0, 1e-5))) {
        #ifndef CUDA
        fprintf(stderr,"m2<0 (%e)\n", g->m2);
        #endif
        *status=0; 
        return;
    }

    g->mK = 1./sqrt(a2*(g->m2p+g->m2m));

//    if (a==0.0){
//        g->m2 = 0.0;
//        g->mK = 1./sqrt(fabs(l2-q));
//    }

    // periastron
    g->rp = (g->nrr > 0) ? creal(g->r1) : 0.0;
    g->rp_R = geodesic_s2i_int_R_periastron(g);
}



DEVICEFUNC
double geodesic_s2i_theta_infinity(double spin, double l, double q, double r, double m)
// returns inclination angle (cosine of theta) at which photon arives to infinity
// inputs are BH spin, motion constants [l,q] and a position at photon trajectory [r,m]
// [r,m] represent an arbitrary position along the photon trajectory 
{
/*
    int status;
    double R, Tpp, Tep, cos_inc;
    geodesic g;
    geodesic_s2i_setup_local(1e-8, p.l, p.q, &g, &status);
    if (!status) return NAN;
    if ((fabs(p.lq_r) < g.rp) && (fabs(p.lq_r)+1e-6>g.rp)) p.lq_r+=sign(p.lq_r)*1e-6;
                if (fabs(p.lq_r) < g.rp) break;//{fprintf(stderr,"XX>lq_r=%.8f  rp=%.8f\n",p.lq_r,g.rp);break;}
                R = geodesic_s2i_int_R(&g, fabs(p.lq_r), (p.lq_r<0));
                Tpp = 2.*g.mK*elliptic_k(g.m2);  // int_{-\mu_plus}^{\mu_plus}
                Tep = g.mK*jacobi_icn(p.lq_m/sqrt(g.m2p),g.m2);     // int_{\mu_em}^{\mu_plus}
                while (R > Tpp) {R -= Tpp;}
                if (R > Tep) {
                    cos_inc = sqrt(g.m2p)*jacobi_cn((R-Tep)/g.mK,g.m2);
                    p.k[2] = +fabs(p.k[2]);
                } else {
                    cos_inc = sqrt(g.m2p)*jacobi_cn((Tep-R)/g.mK,g.m2);
                    p.k[2] = -fabs(p.k[2]);
                }
                p.pos[2] = cos_inc;
*/                
return 0.0;
}

DEVICEFUNC
double geodesic_s2i_int_R_periastron(geodesic *g)
// the value of R-integral at periastron (radial turning point)
{
    double r1,r2,r3,r4,u,v;
    
    if (g->nrr == 4) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        r3 = creal(g->r3);
        r4 = creal(g->r4);

        if (r1==r2) {
            #ifndef CUDA
            fprintf(stderr, "ERR: no solution implemented for r1==r2\n");
            #endif
            return 0.0;
        } else {
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            return 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r2-r4)/(r1-r4)), m4);
        }
    } else
    if (g->nrr == 2) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        u  = creal(g->r3);
        v  = cimag(g->r3);
        double A = sqrt(sqr(r1-u)+sqr(v));
        double B = sqrt(sqr(r2-u)+sqr(v));
        double m2 = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
        return 1./sqrt(A*B) * jacobi_icn((A-B)/(A+B), m2);
    } else {
        #ifndef CUDA
        //fprintf(stderr, "ERR: no solution implemented for r1-r4 all complex\n");
        #endif
        return 0.0;
    }
}



DEVICEFUNC
double geodesic_s2i_int_R(geodesic *g, double r, int beyond_pa)
// the value of R-integral at radius r
{
    double r1,r2,r3,r4,u,v;
    
    if (r  < g->rp) error("geodesic_s2i_int_R: r < periastron (%.3e < %.3e)", r, g->rp);
    if (r == g->rp) return g->rp_R;
    
    if (g->nrr == 4) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        r3 = creal(g->r3);
        r4 = creal(g->r4);

        if (r1==r2) {
            #ifndef CUDA
            fprintf(stderr, "WRN: r1 == r2\n");
            #endif
            return 0.0;
        } else {
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            double R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt(((r2-r4)*(r-r1))/((r1-r4)*(r-r2))), m4);
            return (beyond_pa) ? g->rp_R + R : g->rp_R - R;
        }
    } else
    if (g->nrr == 2) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        u  = creal(g->r3);
        v  = cimag(g->r3);
        double A = sqrt(sqr(r1-u)+sqr(v));
        double B = sqrt(sqr(r2-u)+sqr(v));
        double m2 = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
        double R  = 1./sqrt(A*B) * jacobi_icn(((A-B)*r+r1*B-r2*A)/((A+B)*r-r1*B-r2*A), m2);
        return (beyond_pa) ? g->rp_R + R : g->rp_R - R;
        //return g->rp_R - R;
        //TODO??? nema tu by take return (beyond_pa) ? g->rp_R + R : g->rp_R - R; jako nahore?
    } else {
        #ifndef CUDA
        fprintf(stderr, "ERR: no roots\n");
        #endif
        return 0.0;
    }
}



DEVICEFUNC
double geodesic_s2i_inv_R(geodesic *g, double R, int* beyond_pa)
// radius at which the R-integral gains the value R
{
    double r1,r2,r3,r4,u,v;
    
    if ((R<=0.0)||(R>=2.*g->rp_R)) return -1.0;
    if (R == g->rp_R) return g->rp;
    if (beyond_pa) *beyond_pa = (R > g->rp_R);
    
    if (g->nrr == 4) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        r3 = creal(g->r3);
        r4 = creal(g->r4);

        if (r1==r2) {
            #ifndef CUDA
            fprintf(stderr, "WRN: r1 == r2\n");
            #endif
            return 0.0;
        } else {
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            double x4 = 0.5*(R - g->rp_R)*sqrt((r2-r4)*(r1-r3));
            double sn2 = pow( jacobi_sn(x4,m4), 2.0);
            return ( r1*(r2-r4)-r2*(r1-r4)*sn2 ) / ( r2-r4-(r1-r4)*sn2 );
        }
    } else
    if (g->nrr == 2) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        u  = creal(g->r3);
        v  = cimag(g->r3);
        double A = sqrt(sqr(r1-u)+sqr(v));
        double B = sqrt(sqr(r2-u)+sqr(v));
        double m2 = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
        if ((creal(g->r1)<r_bh(g->a)) && (R > g->rp_R)) {/*fprintf(stderr, "WRN: geodesic_s2i_inv_R - unphysical solution\n");*/ return -1.0;}  //unphysical solution
        double cn = jacobi_cn(sqrt(A*B)*(R - g->rp_R), m2);
        return (r2*A - r1*B - (r2*A+r1*B)*cn ) / ( (A-B) - (A+B)*cn );
    } else {
        #ifndef CUDA
        fprintf(stderr, "ERR: no roots\n");
        #endif
        return 0.0;
    }
}



DEVICEFUNC
double geodesic_s2i_int_T_eqplane(geodesic *g, int order, double *dk2dx)
// the value of T-integral in the equatorial plane
{
    if (g->q<=0.0) {
        #ifndef CUDA
        fprintf(stderr,"WRN: q<0\n"); 
        #endif
        return -1.0; 
    }

    double u = g->cos_i/sqrt(g->m2p);
    if (!ensure_range(&u, 0.0, 1.0, 1e-4)) {
        #ifndef CUDA
        fprintf(stderr,"u out of range (%e)\n", u); 
        #endif
        return -1.0;
    }
    
    if (dk2dx) (*dk2dx) = (order%2==0)?-1.:+1.; 

    if (g->beta == 0.0) return g->mK*elliptic_k(g->m2);

    double sgn_psi = (g->beta>0.0) ? 1.0 : -1.0;
    double psi = jacobi_icn(u,g->m2);

    if (isnan(psi)) {
        #ifndef CUDA
        fprintf(stderr,"WRN: psi is nan - icn(%.7e,%.7e) a=%.2e, b=%.2e, ci=%.2e,m2p=%.2e\n",u, g->m2,g->alpha,g->beta,g->cos_i,sqrt(g->m2p)); 
        #endif
        return -10.0; //psi=0.0; 
    }

    double res = g->mK*( (2.*(double)order+1.)*elliptic_k(g->m2) + sgn_psi*psi );
    //fprintf(stderr,"intT: %d K=%.2e psi=%.2e u=%.2e res=%.2e rpr=%.2e/%.2e\n", order, g->mK*elliptic_k(g->m2), sgn_psi*g->mK*psi, u, res, g->rp_R, g->rp);

    if (isnan(res)) {
        #ifndef CUDA
        fprintf(stderr,"ERR: int_theta is nan (psi=%e, g.m2=%e)\n", psi, g->m2); 
        #endif
        return -1.0; 
    }

    return res;
}



DEVICEFUNC
double geodesic_s2i_inv_T(geodesic *g, double T, double *dk2dx)
// theta angle at which the T-integral gains the value T
// optional output dk2dx defines the sign (+1/-1) of the P^theta component of photon 4-momentum
{
    //if (g->beta == 0.0) return acos(sqrt(g->m2p));

    double Tmo = g->mK*jacobi_icn(g->cos_i/sqrt(g->m2p),g->m2);     // int_{\mu_obs}^{\mu_plus}
    double Tmp = g->mK*elliptic_k(g->m2);                            // int_{0}^{\mu_plus}
    int k = (g->beta>=0.0) ? 3 : 0;
    T += (g->beta>=0.0) ? Tmp-Tmo : Tmo;                             // shift T to the nearest higher \mu=0 or \mu=\mu_plus

    for (;T>Tmp; k++) T -= Tmp;
    switch(k%4) {
        case 0: if (dk2dx) (*dk2dx)=-1; return +sqrt(g->m2p)*jacobi_cn(T/g->mK,g->m2);
        case 1: if (dk2dx) (*dk2dx)=-1; return -sqrt(g->m2p)*jacobi_cn((Tmp-T)/g->mK,g->m2);
        case 2: if (dk2dx) (*dk2dx)=+1; return -sqrt(g->m2p)*jacobi_cn(T/g->mK,g->m2);
        case 3: if (dk2dx) (*dk2dx)=+1; return +sqrt(g->m2p)*jacobi_cn((Tmp-T)/g->mK,g->m2);
    }
    error("geodesic_s2i_inv_T: ???");
    return 0.0;
}


DEVICEFUNC
double geodesic_s2i_timedelay(geodesic *g, double x, double *opt_r, double *opt_m)
{
    int    bpr  = ((g->nrr==4) && (x > g->rp_R));              // beyond periastron radius
    double a2 = sqr(g->a);
    double rp = 1. + sqrt(1.-a2);
    double rm = 1. - sqrt(1.-a2);
    double T  = 0.0;
    double r    = (opt_r) ? *opt_r : geodesic_s2i_inv_R(g,x,&bpr);
    //unused double mu_e = (opt_m) ? *opt_m : geodesic_s2i_inv_T(g,x,NULL);

    double RMAX = 1000.;
    if (r > RMAX) error("geodesic_s2i_timedelay: r > RMAX");
    if (r < g->rp) error("geodesic_s2i_timedelay: r < r_p");

    if (g->nrr == 4) {
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        double r3 = creal(g->r3);
        double r4 = creal(g->r4);
        double R0 = integral_R_r0_re(r1, r2, r3, r4, RMAX) + (bpr?+1:-1)*integral_R_r0_re(r1, r2, r3, r4, r);
        double R1 = integral_R_r1_re(r1, r2, r3, r4, RMAX) + (bpr?+1:-1)*integral_R_r1_re(r1, r2, r3, r4, r);
        double R2 = integral_R_r2_re(r1, r2, r3, r4, RMAX) + (bpr?+1:-1)*integral_R_r2_re(r1, r2, r3, r4, r);
        double RA = integral_R_rp_re(r1, r2, r3, r4, rp, RMAX) + (bpr?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rp, r);
        double RB = integral_R_rp_re(r1, r2, r3, r4, rm, RMAX) + (bpr?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rm, r);
        double a = +(g->a*g->l-4.)*rp + 2.*a2;
        double b = -(g->a*g->l-4.)*rm + 2.*a2;
        T += 4.*R0 + 2.*R1 + R2 - (a*RA + b*RB)/sqrt(1.-a2);
    } else 
    if (g->nrr == 2) {
        //if (bpr) error("geodesic_s2i_timedelay: cannot be (beyond_pa) && (nrr==2)");
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        double R0 = integral_R_r0_cc(r1, r2, g->r3, RMAX) - integral_R_r0_cc(r1, r2, g->r3, r);
        double R1 = integral_R_r1_cc(r1, r2, g->r3, r, RMAX);// - integral_R_r1_cc(r1, r2, g->r3, r);
        double R2 = integral_R_r2_cc(r1, r2, g->r3, r, RMAX);// - integral_R_r2_cc(r1, r2, g->r3, r);
        double RA = integral_R_rp_cc2(r1, r2, g->r3, rp, r, RMAX);// + integral_R_rp_cc(r1, r2, g->r3, rp, r);
        double RB = integral_R_rp_cc2(r1, r2, g->r3, rm, r, RMAX);// + integral_R_rp_cc(r1, r2, g->r3, rm, r);
        double a = +(g->a*g->l-4.)*rp + 2.*a2;
        double b = -(g->a*g->l-4.)*rm + 2.*a2;
        T += 4.*R0 + 2.*R1 + R2 - (a*RA + b*RB)/sqrt(1.-a2);
    } else {
        // not implemented
    }
    
/* 
    double Tmp = g->mK*elliptic_k(g->m2);
    double Tmo = g->mK*jacobi_icn(g->cos_i/sqrt(g->m2p),g->m2);
    double Pmp = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, 0.0);  
    double Pmo = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, g->cos_i); 
    double Pme = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, fabs(mu_e));

    x += (g->beta>=0.0) ? Tmp-Tmo : Tmo;                             // shift T to the nearest higher \mu=0 or \mu=\mu_plus
    T -= (g->beta>=0.0) ? Pmp-Pmo : Pmo;
    int k = (g->beta>=0.0) ? 3 : 0;
    while (x >= Tmp) {
        x -= Tmp;
        T += Pmp;
        k++;
    }
    switch(k%4) {
        case 0: T += Pme; break;
        case 1: T += Pmp-Pme; break;
        case 2: T += Pme; break;
        case 3: T += Pmp-Pme; break;
    }
 */
    return T-RMAX;
}


DEVICEFUNC
double geodesic_s2i_phi(geodesic *g, double x, double *opt_r, double *opt_m)
{
    int    bpr  = (g->nrr==4) && (x > g->rp_R);              // beyond periastron radius
    double a2 = sqr(g->a);
    double rp   = 1. + sqrt(1.-a2);
    double rm   = 1. - sqrt(1.-a2);
    double phi  = 0.0;
    double r    = (opt_r) ? *opt_r : geodesic_s2i_inv_R(g,x,&bpr);
    double mu_e = (opt_m) ? *opt_m : geodesic_s2i_inv_T(g,x,NULL);

    // part 1 - R integral
    if (g->nrr == 4) {
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        double r3 = creal(g->r3);
        double r4 = creal(g->r4);
        double A = integral_R_rp_re_inf(r1, r2, r3, r4, rp) + (bpr?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rp, r); 
        double B = integral_R_rp_re_inf(r1, r2, r3, r4, rm) + (bpr?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rm, r);
        phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
    } else 
    if (g->nrr == 2) {
        if (bpr) error("geodesic_s2i_phi: cannot be (beyond_pa) && (nrr==2)");
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        double A = integral_R_rp_cc2_inf(r1, r2, g->r3, rp, r); 
        double B = integral_R_rp_cc2_inf(r1, r2, g->r3, rm, r); 
        phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
    } else {
        error("geodesic_s2i_phi: g->nrr != [2,4]");
    }

    double Tmp = g->mK*elliptic_k(g->m2);
    double Tmo = g->mK*jacobi_icn(g->cos_i/sqrt(g->m2p),g->m2);
    double Pmp = -g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, 0.0);         // the front minus should be understood
    double Pmo = -g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, g->cos_i);   //  -- " --
    double Pme = -g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, fabs(mu_e));  //  -- " --

    x   += (g->beta>=0.0) ? Tmp-Tmo : Tmo;                             // shift T to the nearest higher \mu=0 or \mu=\mu_plus
    phi -= (g->beta>=0.0) ? Pmp-Pmo : Pmo;
    int k = (g->beta>=0.0) ? 3 : 0;
    while (x >= Tmp) {
        x -= Tmp;
        phi += Pmp;
        k++;
    }
    switch(k%4) {
        case 0: phi += Pme; break;
        case 1: phi += Pmp-Pme; break;
        case 2: phi += Pme; break;
        case 3: phi += Pmp-Pme; break;
    }

    while (phi > 2.*M_PI) phi -= 2.*M_PI;
    while (phi < 0.0) phi += 2.*M_PI;
    return phi;
}


DEVICEFUNC
double geodesic_s2i_afp(geodesic *g, double x, double *opt_r, double *opt_m)
// affine parameter
{
    int    bpr  = (g->nrr==4) && (x > g->rp_R);              // beyond periastron radius
    double afp  = 0.0;
    double r    = (opt_r) ? *opt_r : geodesic_s2i_inv_R(g,x,&bpr);
    double mu_e = (opt_m) ? *opt_m : geodesic_s2i_inv_T(g,x,NULL);

    double RMAX = 1000.;
    if (g->nrr == 4) {
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        double r3 = creal(g->r3);
        double r4 = creal(g->r4);
        afp += integral_R_r2_re(r1, r2, r3, r4, RMAX) + (bpr?+1:-1)*integral_R_r2_re(r1, r2, r3, r4, r);
    } else 
    if (g->nrr == 2) {
        //if (bpr) error("geodesic_s2i_timedelay: cannot be (beyond_pa) && (nrr==2)");
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        afp += integral_R_r2_cc(r1, r2, g->r3, r, RMAX);
    } else {
        error("geodesic_s2i_afp: g->nrr != [2,4]");
    }

    double Tmp = g->mK*elliptic_k(g->m2);
    double Tmo = g->mK*jacobi_icn(g->cos_i/sqrt(g->m2p),g->m2);
    double Pmp = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, 0.0);  
    double Pmo = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, g->cos_i); 
    double Pme = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, fabs(mu_e));

    x   += (g->beta>=0.0) ? Tmp-Tmo : Tmo;                             // shift T to the nearest higher \mu=0 or \mu=\mu_plus
    afp -= (g->beta>=0.0) ? Pmp-Pmo : Pmo;
    int k = (g->beta>=0.0) ? 3 : 0;
    while (x >= Tmp) {
        x -= Tmp;
        afp += Pmp;
        k++;
    }
    switch(k%4) {
        case 0: afp += Pme; break;
        case 1: afp += Pmp-Pme; break;
        case 2: afp += Pme; break;
        case 3: afp += Pmp-Pme; break;
    }

    return afp-RMAX;
}


DEVICEFUNC
double geodesic_s2i_mue_eqplane(double k[4], double U[4], double n[4], 
       double gtt, double grr, double gmm, double gff, double gtf)
// cosine of the emmision angle
// mue = cos(delta) = -(p.n)/(p.U)
{
    double A = gtt*k[0]*n[0]+grr*k[1]*n[1]+gmm*k[2]*n[2]+gff*k[3]*n[3]+gtf*k[0]*n[3]+gtf*k[3]*n[0];
    double B = gtt*k[0]*U[0]+grr*k[1]*U[1]+gmm*k[2]*U[2]+gff*k[3]*U[3]+gtf*k[0]*U[3]+gtf*k[3]*U[0];
    return -A/B;
}


DEVICEFUNC
void geodesic_s2i_solution_eqplane(geodesic *g, int order, double *r, int *beyond_pa, double *phi, double *time_delay, int *status)
{
    double ksgn_2;
    (*status) = 0;
    if (g->q < 0.0) return;
    g->x = geodesic_s2i_int_T_eqplane(g, order, &ksgn_2);
    (*r) = geodesic_s2i_inv_R(g, g->x, beyond_pa);
    if ((*r) < r_bh(g->a)) return;    
    (*status) = 1;    
    photon_momentum(g->a, *r, 0.0, g->l, g->q, (*beyond_pa)?-1:+1, ksgn_2, g->k);

    double m = 0.0;
    if (time_delay) *time_delay = geodesic_s2i_timedelay(g, g->x, r, &m);
    if (phi) *phi = geodesic_s2i_phi(g, g->x, r, &m);
//if (beyond_pa) *status=0;
}



/*
void geodesic_s2i_solution_surface(geodesic *g, double *r, double *z, double (*H)(double), int *status)
{
    double rh = r_bh(g->a);

    double fx(geodesic *g, double T) {
        int pax;
        double _r = geodesic_s2i_inv_R(g, T, &pax); if (_r<0.0) {fprintf(stderr,"err: inv_R<0 (T=%e _r=%e Trp=%e)\n", T, _r, g->rp_R); exit(-1);}
        double _t = geodesic_s2i_inv_T(g, T, NULL);
        double _z = (_r*cos(_t));
        double _x = sqrt(_r*_r - _z*_z);
        double _H = (*H)(_x);
        if (_H/_x < 1e-5) _H = 1e-5*_x;
        //double surface_z = 0.0;//(_r>6.0) ? 15.*atan((_r-6.0)/8.)/3.1415*2.0 : 1.0;
        //fprintf(stderr,"r=%e  z=%e  H=%e\n", _r, _z, (*H)(_r));
        return (_z - _H);
        
    }

    long gdbis(double x1, double x2, double xacc, double (*fx)(geodesic*,double), double* result)
    {
        double dx, f, fmid, xmid, rtb;
        long j=0;
        int MAX_STEPS = 500;

        // what is it doing?    
        if (geodesic_s2i_inv_R(g,x2,NULL)<0.0){
            #ifndef CUDA
            fprintf(stderr,"WRN: x2 adjustment in gdbis: (x1=%e,r=%e) x2=%e --> ", x1,geodesic_s2i_inv_R(g,x2,NULL),x2);
            //x2 = geodesic_s2i_int_R(g,g->r1,0);
            fprintf(stderr,"%e\n", x2);
            #endif
        }

        fmid = (*fx)(g,x2);
        f    = (*fx)(g,x1);
        if ((f*fmid) >= 0.0) {
            //fprintf(stderr,"rtbis: root is not bracketed %e/%e", f,fmid);
            //double _r;
            //fprintf(stderr,"z1=%e, x2=%e, r=%e\n",fmid,x2,geodesic_s2i_inv_R(g,x2,NULL));
            //fprintf(stderr,"z2=%e, x1=%e, r=%e\n",f,x1,geodesic_s2i_inv_R(g,x1,NULL));
            return 0;
        } //{fprintf(stderr,"nob\n");return(0);}//error("rtbis: root is not bracketed");
    
        if (f < 0.0) {
            rtb = x1;
            dx  = x2-x1;
        }
        else {
            rtb = x2;
            dx  = x1-x2;
        }

        for (j=0; j<MAX_STEPS; j++) {
            dx = dx*0.5;
            xmid = rtb+dx;
            fmid = (*fx)(g,xmid);
            if (fmid <= 0.0) rtb = xmid;
            if ((fabs(dx) < xacc) || (fmid = 0.0)) break;
        }
        if (j >= MAX_STEPS) error("rtbis: too many steps"); 
    
        *result = rtb;
        return(1);
    }
    //---

    *status = 0;
    double ksgn_2;
    double T0;
    double T1;
    int pax;
    double result;

    T0 = g->rp_R>1e-3 ? 1e-4 : g->rp_R/10.;
    T1 = geodesic_s2i_int_T_eqplane(g,0,&ksgn_2);
//    if (g->nrr >= 4) {

    *r = geodesic_s2i_inv_R(g, T1, &pax);
    if ((T1<0) || (*r < rh)) T1 = geodesic_s2i_int_R(g, rh*1.01, 0);
    
    if ((*H)(*r) == 0.0) {
        *z = 0.0;
        g->x = T1;
        photon_momentum(g->a, g->q, g->l, *r, 0.0, (T1>g->rp_R)?-1:+1, ksgn_2, g->k);
        *status = 1;
        return;
    } 

    // check if the outer anchor point is above the surface
    *r = geodesic_s2i_inv_R(g, T0, &pax);
    *z = (*r)*geodesic_s2i_inv_T(g, T0, &ksgn_2);
    if ((*H)(sqrt(sqr(*r)-sqr(*z))) >= *z) {
        //fprintf(stderr, "geodesic_s2i_solution_surface:: outer anchor inside disk (r=%e z=%e, H=%e)\n", *r, *z, (*H)(sqrt(sqr(*r)-sqr(*z))));
        *status = 0;
        return;
    } 

    *status = (int)gdbis(T0, T1, 1e-8, &fx, &result);
    if (!(*status)) return;

    if ((result<T0) || (result>T1)) {
        #ifndef CUDA
        fprintf(stderr, "yygdbis(T0=%e, T1=%em res=%e)\n", T0, T1, result);
        #endif
        *status = 0;
        return;
    } 

    if (result == T1) {
        *z = 0.0;
        g->x = result;
        photon_momentum(g->a, g->q, g->l, *r, 0.0, (T1>g->rp_R)?-1:+1, ksgn_2, g->k);
        *status = 1;
        return;
    } 

    *r = geodesic_s2i_inv_R(g, result, &pax);

    if ((!(*status)) || (*r<rh)) {
        *status = 1;
        return;
    } 

//    } // if q>0
//    else {
//      result = geodesic_s2i_int_R(g, rh*1.05, 0);
//      if (result<0.0)fprintf(stderr, "q>0 %e q=%e\n", result, g->q);
//    }

    if (acos(g->cos_i) > deg2rad(45.)&&(T0>9e-5)) {
        double _r, _z, _R, _x, _h;
        _x = 0.998*result;
        do {
            _r = geodesic_s2i_inv_R(g, _x, &pax);
            _z = (_r)*geodesic_s2i_inv_T(g, _x, &ksgn_2);
            _R = sqrt(_r*_r - _z*_z);
            _h = (*H)(_R);
            if ((_z < _h)) {
                //fprintf(stderr, "xxgdbis(T0=%e, T1=%em r=%e, z=%e, h=%e)\n", _x/1e4, _x, _r, _z, _h);
                //if (g->q < 0) fprintf(stderr, "adjust surface: x=%e  r=%e z=%e h=%e --> ", _x, _r, _z, _h);
                //fprintf(stderr,"x=%e/%e/%e, r=%e z=%e R=%e h=%e\n", _x,result,g->rp_R, _r, _z, _R, _h);
                *status = (int)gdbis(_x/1e4, _x, 1e-8, &fx, &result);
                if (!(*status)) return;
                _x = result;
                _r = geodesic_s2i_inv_R(g, _x, &pax);
                _z = (_r)*geodesic_s2i_inv_T(g, _x, NULL);
                //if (g->q < 0) fprintf(stderr, "r,z=%e,%e\n", _r, _z);
            }
            _x *= 0.990;
            if (_r>1e4) break;
        }
        while ((_r<4.*r_ms(g->a)) || (_z < 2.*_h)) ;
    }
    
    *r = geodesic_s2i_inv_R(g, result, &pax);
    //if (*r<g->rh) fprintf(stderr,"err: r < rh %e (%e : %e - %e) %e/%e\n",*r,result,geodesic_s2i_inv_R_dbg(g,T1,NULL),geodesic_s2i_inv_R_dbg(g,T0,NULL),g->rp_R,g->rp_T);
    *z = (*r)*geodesic_s2i_inv_T(g, result, &ksgn_2);
    g->x = result;
    photon_momentum(g->a, g->q, g->l, *r, (*z)/(*r), ((g->x)>(g->rp_R))?-1:+1, ksgn_2, g->k);
    *status = (*status) && (*r>rh);
}
*/


DEVICEFUNC
void geodesic_s2i_solution_surface(geodesic *g, double *r, double *m, double (*H)(double), double accuracy, int *status)
{
    *r = 0.0; *m = 0.0;
    double rh = r_bh(g->a);
    double r0 = max(1000.0, g->rp*2.);
    double r1, m1, x1, R1, h1, H1;
    double r2, m2, x2, R2, h2, H2;

    do {
      // initialize photon-follow from r0
      geodesic_s2i_follow_init(g, r0, status);
      if (!(*status)) return;

      // get initial position
      x1 = x2 = g->x;
      r1 = r2 = geodesic_s2i_inv_R(g, g->x, NULL);
      m1 = m2 = geodesic_s2i_inv_T(g, g->x, NULL);
      R1 = R2 = r1*sqrt(1.-m1*m1);
      h1 = h2 = r1*m1;
      H1 = H2 = fabs((*H)(R1));

      r0 *= 2.0;
    } while ((h1<=H1) && (r0<5e4));

    // test if the initial point is already under surface
    if (h1 <= H1) { *status=0; return; }

    // follow trajectory with adaptive step until surface is reached
    do {
        double step;
        r2=r1; m2=m1; R2=R1; x2=x1; h2=h1; H2=H1;

        // new step
        step=max((h1-H1)/5., H1*accuracy);
        geodesic_s2i_follow(g, step, &r1, &m1, NULL, NULL, status);
        //fprintf(stderr,"r=%.3f  step=%.3e H1=%.3e, sta=%d\n", r1, step, H1, *status);
    	if (!*status) return;

        // current photon and surface location
        x1 = g->x;
        R1 = r1*sqrt(1.-m1*m1);
        h1 = r1*m1;
        H1 = fabs((*H)(R1));

        // surface hit?
        if (h1 <= H1) {
            *r = r2;
            *m = m2;
            *status = 1;
            return;
        }

        // eq plane hit?
        if (fabs(h1) < accuracy) {
            *r = r1;
            *m = m1;
            *status = 1;
            return;
        }

        // terminating conditions
        //if (r1>1.1*r0) break;
        if (m1<0.0) break;
        if (r1<1.05*rh) break;
	} while (*status);

    *status = 0;
}


DEVICEFUNC
void geodesic_s2i_follow_init(geodesic *g, double rmax, int *status)
{
/*
    #ifndef CUDA
    if (rmax<10.0) fprintf(stderr, "WRN: rmax should be >10 (geodesic_s2i_follow_init)\n");
    #endif
    double x1,x2,t1,t2;
    *status = 0;
    if (rmax < g->rp) return;
    x1 = geodesic_s2i_int_R(g, rmax, 0);
    x2 = geodesic_s2i_int_R(g, rmax*1.01, 0);
    t1 = geodesic_s2i_timedelay(g, x1, NULL, NULL);
    t2 = geodesic_s2i_timedelay(g, x2, NULL, NULL);
    g->x     = x1;
    g->_t    = t1;
    g->_dxds = (x2-x1)/(t2-t1);
    *status = 1;
//    fprintf(stderr, "I> x=%.3e  t=%.3e  dxds=%.3e\n", g->x, g->_t, g->_dxds);
//    fprintf(stderr, "I> f=%.3f\n", geodesic_s2i_phi(g, g->x, NULL, NULL)*180/M_PI);
*/
    #ifndef CUDA
    if (rmax<10.0) fprintf(stderr, "WRN: rmax should be >10 (geodesic_s2i_follow_init)\n");
    #endif
    *status = 0;
    if (rmax < g->rp) return;
    g->x = geodesic_s2i_int_R(g, rmax, 0);
    *status = 1;
}


DEVICEFUNC
void geodesic_s2i_follow(geodesic *g, double step, double *r, double *m, double *f, double *t, int *status)
{
/*
    int brp;
    double ksgn_2;
    double dx = step*g->_dxds;

    //do_step_fwd:
    double new_x = g->x+dx;
//    fprintf(stderr, "F> x=%.5e  newx=%.5e  dx=%.5e\n", g->x, new_x, dx);
    double new_r = geodesic_s2i_inv_R(g, new_x, &brp);
    if (new_r < 1.05*r_bh(g->a)) {
        *status = 0;
        return;
    }
//    fprintf(stderr, "F> new_r=%.5f (brp=%d)\n", new_r, brp);
    double new_t = geodesic_s2i_timedelay(g, new_x, &new_r, NULL);
//    fprintf(stderr, "F> new_t=%.5f\n", new_t);

    double new_ds = fabs(new_t - g->_t);
    dx  = dx * (step/new_ds);
//    if (fabs(1.-(ds/step))>0.05) {fprintf(stderr,"WRN: step difference > 5%%\n");/ * goto do_step_fwd;* /}

    //if (dx<0.0) {fprintf(stderr, "WRN: dx<0.0 (t1=%e  t2=%e  dx=%e  r=%e)\n", dt, dt2, dx, r);break;}

    g->x     = new_x;
    g->_t    = new_t;
    g->_dxds = dx/new_ds;

    if (r) (*r) = new_r;
    if (m) (*m) = geodesic_s2i_inv_T(g, new_x, &ksgn_2);
//    fprintf(stderr, "F> m=%.5f\n", acos(*m)*180./M_PI);
    if (f) (*f) = geodesic_s2i_phi(g, new_x, r, m);
//    fprintf(stderr, "F> r=%.2f  f=%.5f (%d)\n", *r, *f*180./M_PI, brp);
    if (t) (*t) = new_t;
    if (ds)(*ds)= new_ds;
    if ((r)&&(m)) photon_momentum(g->a, g->q, g->l, *r, *m, brp?-1:+1, ksgn_2, g->k);
    *status = (new_r > r_bh(g->a));
*/

    int brp;
    double ksgn_2;
    const double MAXSTEP = 10.0;

    double _x = g->x;
    double _r = (*r);
    double _m = (*m);
    
    do {
        double truestep = min(step, MAXSTEP);
        _x = _x + truestep/(sqr(_r)+sqr(g->a*_m));   // d(afp)/d(x) = r^2 + a^2*m^2
        _r = geodesic_s2i_inv_R(g, _x, &brp);
        _m = geodesic_s2i_inv_T(g, _x, &ksgn_2);
        if (_r < 1.01*r_bh(g->a)) {
            *status = 0;
            return;
        }
        step -= truestep;
    } while (step > 1e-5);

    g->x = _x;
    photon_momentum(g->a, _r, _m, g->l, g->q, brp?-1:+1, ksgn_2, g->k);
    if (r) (*r) = _r;
    if (m) (*m) = _m;
    if (f) (*f) = geodesic_s2i_phi(g, _x, &_r, &_m);
    if (t) (*t) = geodesic_s2i_timedelay(g, _x, &_r, &_m);

    *status = 1;
}


