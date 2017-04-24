//************************************************************************
//    SIM5 library
//    sim5kerr-geod.c - null-geodesic motion routines
//------------------------------------------------------------------------
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (c) Michal Bursa, Astronomical Institute of the CAS
//************************************************************************

/*
#ifdef CUDA
__host__ __device__ void error(char *s) {
    //abort();
}
#endif
*/

// define some helper macros
#define theta_int(x) (g->mK*jacobi_icn((x)/sqrt(g->m2p),g->m2))
#define theta_inv(x) (sqrt(g->m2p)*jacobi_cn((x)/g->mK,g->m2))

// unit-private function declarations
DEVICEFUNC int geodesic_priv_R_roots(geodesic *g, int *status);
DEVICEFUNC int geodesic_priv_T_roots(geodesic *g, int *status);



DEVICEFUNC
int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *status)
//! Makes setup for geodesic that is specified by impact parameters at infinity.
//!
//! Parameters:
//!     i      - inclination angle of observer (angle between BH rotation axis and
//!              direction to observer) [radians]
//!     a      - BH spin [0..1]
//!     alpha  - impact parameter in horizontal direction [GM/c^2]
//!     beta   - impact parameter in vertical direction [GM/c^2]
//!     g      - structure with information about geodesic
//!     status - status variable
//!
//! Returns:
//!     * status is 1 if setup is sucessfull, 0 in case of an error
//!     * information about geodesic is stored in structure g
{
    if (a < 1e-8) a = 1e-8;

    g->a = a;
    g->i = i;
    g->cos_i = cos(i);
    g->alpha = alpha;
    g->beta  = beta;
    g->x = 0.0;

    // constants of motion
    g->l = -alpha*sin(i);
    g->q = sqr(beta) + sqr(cos(i))*(sqr(alpha)-sqr(a));

    if (g->q == 0.0) {
        // we are not yet ready to handle this case (Dexter&Agol provide a solution)
        // skip it silently
        if (status) *status = 0;
        #ifndef CUDA
        fprintf(stderr,"q=0\n");
        #endif
        return FALSE;
    }

    // get geodesic
    if (!geodesic_priv_R_roots(g, status)) return FALSE;
    if (!geodesic_priv_T_roots(g, status)) return FALSE;

    // value of T-integral between turning points \int[-\mu_plus..\mu_plus]
	g->Tpp = 2.*theta_int(0.0);

    // value of T-integral between observer's position and turning point \int[cos_i..\mu_plus]
	g->Tip = theta_int(g->cos_i);

    if (status) *status = 1;
    return TRUE;
}




DEVICEFUNC
int geodesic_init_src(double a, double r, double m, double k[4], int bpa, geodesic *g, int *status)
//! Makes setup for geodesic that is specified by a point and direction (4-momentum vector).
//!
//! Parameters:
//!     a      - BH spin [0..1]
//!     r      - radial coordinate [GM/c^2]
//!     m      - poloidal coordinate [cos(theta)]
//!     k      - direction of the photon (4-momentum vector)
//!     bpa    - position with respect to periastron (0=before periastron, 1=after periastron)
//!     g      - structure with information about geodesic
//!     status - status variable
//!
//! Returns:
//!     * status is 1 if setup is sucessfull, 0 in case of an error
//!     * information about geodesic is stored in structure g
{
    if (a < 1e-8) a = 1e-8;

    // calculate motion constants
    double l,q;
    photon_motion_constants(a, r, m, k, &l, &q);

    g->a = a;
    g->l = l;
    g->q = q;
    g->x = 0.0;

    // ...
    g->i = g->cos_i = g->alpha = g->beta = NAN;

    if (!geodesic_priv_R_roots(g, status)) return FALSE;
    if (!geodesic_priv_T_roots(g, status)) return FALSE;
    
    // r is bellow radial turning point - this geodesic cannot escape
    if (r < g->rp) {
        if (status) *status = 0;
        return FALSE;
    }

    // determine at-infinity parameters of the geodesic (inclination, impact parameters)
    // - inclination is poloidal coordinate to which the photon in its current direction will arrive
    if (isnan(g->cos_i)) {
        double T, Tmp, Tpp, sign_dm;
        Tmp = theta_int(m);                         // int_{\mu}^{\mu_plus}
        Tpp = 2.*theta_int(0.0);                    // int_{-\mu_plus}^{\mu_plus}
        T = geodesic_P_int(g, r, bpa);              // int_{r}^{\inf}
        sign_dm = (k[2]<0.0) ? +1.0 : -1.0;           // increasing m = decreasing theta
        //fprintf(stderr,"Tmp=%.3e, Tpp=%.3e, T=%.3e, sdm=%+.0e\n", Tmp, Tpp, T, sign_dm);
        T += (sign_dm>0.0) ? Tpp-Tmp : Tmp;
        //fprintf(stderr,"T2=%.3e\n", T);
        while (T > Tpp) {
            T -= Tpp;
            sign_dm = -sign_dm;
        }
        //fprintf(stderr,"T3=%.3e\n", T);
        g->cos_i = -sign_dm*theta_inv(T);
        g->i     = acos(g->cos_i);
        g->alpha = -g->l/sqrt(1.0-sqr(g->cos_i));
        g->beta  = -sign_dm * sqrt(g->q - sqr(g->cos_i)*(sqr(g->alpha)-sqr(g->a)));
    }

    // value of T-integral between turning points \int[-\mu_plus..\mu_plus]
	g->Tpp = 2.*theta_int(0.0);

    // value of T-integral between observer's position and turning point \int[cos_i..\mu_plus]
	g->Tip = theta_int(g->cos_i);

    if (status) *status = 1;
    return TRUE;
}




DEVICEFUNC
double geodesic_P_int(geodesic *g, double r, int bpa)
//! Returns the value of position integral along geodesics at point r.
//!
//! The function gives the value of the integral
//!     P = \int 1/\sqrt{R} dr = \int 1/\sqrt{\Theta} d\theta
//! This integral is integrated from infinity to the point, where the geodesic
//! reaches radius r either before or behind its periastron.
//! The value of the integral increases monotonicly from infinity along the geodesic
//! and thus we use it for parametrizing it.
//! Note it is not affine parameter, which would be other choice for parametrization.
//!
//! Parameters:
//!     g      - geodesic
//!     r      - radial coordinate [GM/c^2]
//!     bpa    - position with respect to periastron (0=before periastron, 1=after periastron)
//!
//! Returns:
//!     Value of the position integral between infinity and given point.
{
    double r1,r2,r3,r4,u,v;

    #ifdef CUDA
    if (r  < g->rp) asm("exit;");
    #else
    if (r  < g->rp) error("geodesic_P_int: r < periastron (%.3e < %.3e; nrr=%d)", r, g->rp, g->nrr);
    #endif
    if (r == g->rp) return g->Rpa;

    if (g->nrr == 4) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        r3 = creal(g->r3);
        r4 = creal(g->r4);

        if (r1==r2) {
            #ifndef CUDA
            error("(geodesic_affine): no solution implemented for r1==r2\n");
            #endif
            return 0.0;
        } else {
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            double R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt(((r2-r4)*(r-r1))/((r1-r4)*(r-r2))), m4);
            return (bpa) ? g->Rpa + R : g->Rpa - R;
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
        return (bpa) ? g->Rpa + R : g->Rpa - R;
    } else {
        #ifndef CUDA
        error("no solution implemented for r1-r4 all complex\n");
        #endif
        return 0.0;
    }
}




DEVICEFUNC
void geodesic_position(geodesic *g, double P, double x[4])
//! Returns point on the geodesic, where the position integral gains value P.
//!
//! The integral is evaluted in the form
//! phi = a \int (2r-al)/(Delta \sqrt{R}) dr  +  l \int sin^-2(\theta) / \sqrt{Theta} d\theta
//!
//!
//! Parameters:
//!     g      - geodesic
//!     P      - value of the position integral
//!     x[out] - coordinate
//!
//! Returns:
//!     Fills x[] with position.
{
    return;
}




DEVICEFUNC
double geodesic_position_rad(geodesic *g, double P)
//! Gives radius at which the position integral gains value P.
//!
//! Parameters:
//!     g      - geodesic
//!     P      - value of the position integral
//!
//! Returns:
//!     Radius [GM/c^2]
{
    double r1,r2,r3,r4,u,v;

    if ((P<=0.0)||(P>=2.*g->Rpa)) {
        #ifndef CUDA
        warning("(geodesic_position_rad) P out of range (%e, 2Rpa=%e)\n", P, 2*g->Rpa);
        #endif
        return NAN;
    }
    if (P == g->Rpa) return g->rp;
    //if (beyond_pa) *beyond_pa = (R > g->Rpa);

    if (g->nrr == 4) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        r3 = creal(g->r3);
        r4 = creal(g->r4);

        if (r1==r2) {
            #ifndef CUDA
            warning("(geodesic_position_rad) r1 == r2\n");
            #endif
            return 0.0;
        } else {
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            double x4 = 0.5*(P - g->Rpa)*sqrt((r2-r4)*(r1-r3));
            double sn2 = pow( jacobi_sn(x4,m4), 2.0);
            return ( r1*(r2-r4)-r2*(r1-r4)*sn2 ) / ( r2-r4-(r1-r4)*sn2 );
        }
    } else
    if (g->nrr == 2) {
        if ((creal(g->r1) < r_bh(g->a)) && (P > g->Rpa)) {
            // unphysical solution - geodesic goes to observer through BH
            return NAN;
        }

        r1 = creal(g->r1);
        r2 = creal(g->r2);
        u  = creal(g->r3);
        v  = cimag(g->r3);
        double A = sqrt(sqr(r1-u)+sqr(v));
        double B = sqrt(sqr(r2-u)+sqr(v));
        double m2 = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
        double cn = jacobi_cn(sqrt(A*B)*(P - g->Rpa), m2);
        return (r2*A - r1*B - (r2*A+r1*B)*cn ) / ( (A-B) - (A+B)*cn );
    } else {
        #ifndef CUDA
        warning("no solution implemented for r1-r4 all complex\n");
        #endif
        return 0.0;
    }
}




DEVICEFUNC
double geodesic_position_pol(geodesic *g, double P)
//! Gives poloidal (theta) angle at which the position integral gains value P.
//!
//! Parameters:
//!     g      - geodesic
//!     P      - value of the position integral
//!
//! Returns:
//!     Cosine of theta angle.
{
    //if (g->beta == 0.0) return acos(sqrt(g->m2p));

    double sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
    double T = (sign_dm>0.0) ? -(g->Tpp-g->Tip) : -(g->Tip);

    while (P > T+g->Tpp) {
        T += g->Tpp;
        sign_dm = -sign_dm;
        //fprintf(stderr,"P=%.3e T=%.3e T+Tpp=%.3e delta=%.3e\n", P, T, T+g->Tpp, P-T-g->Tpp);
    }

    return -sign_dm*theta_inv(P-T);
}




DEVICEFUNC
double geodesic_position_azm(geodesic *g, double r, double m, double P)
//! Gives azimuthal (phi) angle at which the position integral gains value P.
//!
//! The value of azimuthal angle is assumed to be zero at infinity and
//! the function gives the change of the angle between
//! the point [r,m] and infinity.
//!
//! Parameters:
//!     g      - geodesic
//!     P      - value of the position integral
//!
//! Returns:
//!     Phi angle in radians (can be more than 2pi)
{
    double phi = 0.0;

    int    bpa  = (g->nrr>0) && (P > g->Rpa);              // beyond periastron
    double a2 = sqr(g->a);
    double rp   = 1. + sqrt(1.-a2);
    double rm   = 1. - sqrt(1.-a2);

    // part 1 - R integral
    if (g->nrr == 4) {
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        double r3 = creal(g->r3);
        double r4 = creal(g->r4);
        double A = integral_R_rp_re_inf(r1, r2, r3, r4, rp) + (bpa?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rp, r);
        double B = integral_R_rp_re_inf(r1, r2, r3, r4, rm) + (bpa?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rm, r);
        phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
    } else
    if (g->nrr == 2) {
        #ifdef CUDA
        if (bpa) asm("exit;");
        #else
        if (bpa) {
            error("(geodesic_position_azm) cannot be (bpa) && (nrr==2)");
            return NAN;
        }
        #endif
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        double A = integral_R_rp_cc2_inf(r1, r2, g->r3, rp, r);
        double B = integral_R_rp_cc2_inf(r1, r2, g->r3, rm, r);
        phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
    } else {
        #ifdef CUDA
        asm("exit;");
        #else
        error("(geodesic_position_azm) g->nrr != [2,4]");
        #endif
    }

    // part 2 - T integral
    // TODO it should be checked if that works for photons with m2m<0
    double phi_pp = 2.0*g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, 0.0);
    double phi_ip =     g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, g->cos_i);
    double phi_mp =     g->l/g->a*integral_T_mp(g->m2m, g->m2p, 1.0, m);

    double T;
    double sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
    if (sign_dm > 0.0) {
        T = -(g->Tpp-g->Tip);
        phi -= phi_pp-phi_ip;
    } else {
        T = -g->Tip;
        phi -= phi_ip;
    }

    while (P >= T+g->Tpp) {
        T += g->Tpp;
        phi += phi_pp;
        sign_dm = -sign_dm;
        break;
    }

    phi += (sign_dm<0) ? phi_mp : phi_pp-phi_mp;
    //phi += (sign_dm<0) ? phi_mp : -phi_mp; //seems wrong

    // change of sign here converts the integration from [r,m] to infinity
    // to integration from infinity to [r,m]
    return phi;
}



DEVICEFUNC
double geodesic_timedelay(geodesic *g, double P1, double P2)
//! Gives travel-time (timedelay) between positions P1 and P2.
//! 
//! Returned value is always positive (independent of relative position of 
//! P1 and P2 along the geodesic).
//!
//! Parameters:
//!     g      - geodesic
//!     P1     - value of the position integral at point A
//!     P2     - value of the position integral at point B
//!
//! Returns:
//!     timedelay (positive) between position P1 and P2 (point A and B)
{
    #ifndef CUDA
    warning("(geodesic_timedelay): not implemented yet\n");
    #endif
    return 0.0;
/*
    double time = 0.0;

    int    bpa  = (g->nrr>0) && (P > g->Rpa);              // beyond periastron
    double a2 = sqr(g->a);
    double rp   = 1. + sqrt(1.-a2);
    double rm   = 1. - sqrt(1.-a2);

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
        time += 4.*R0 + 2.*R1 + R2 - (a*RA + b*RB)/sqrt(1.-a2);
    } else
    if (g->nrr == 2) {
        fprintf(stderr, "ERR (geodesic_position_time): no solution implemented for nrr=2\n");
        / *
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
        time += 4.*R0 + 2.*R1 + R2 - (a*RA + b*RB)/sqrt(1.-a2);
        * /
    } else {
        fprintf(stderr, "ERR (geodesic_position_time): no solution implemented for nrr=0\n");
        // not implemented
    }

    // part 2 - T integral
    // TODO it should be checked if that works for photons with m2m<0
    double time_pp = 2.0*g->l/g->a*integral_T_m2(g->m2m, g->m2p, 1.0, 0.0);
    double time_ip =     g->l/g->a*integral_T_m2(g->m2m, g->m2p, 1.0, g->cos_i);
    double time_mp =     g->l/g->a*integral_T_m2(g->m2m, g->m2p, 1.0, m);

    double T;
    double sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
    if (sign_dm > 0.0) {
        T = -(g->Tpp-g->Tip);
        time -= time_pp-time_ip;
    } else {
        T = -g->Tip;
        time -= time_ip;
    }

    while (P >= T+g->Tpp) {
        T += g->Tpp;
        time += time_pp;
        sign_dm = -sign_dm;
        break;
    }

    time += (sign_dm<0) ? time_mp : time_pp-time_mp;

/ *
    double Tmp = g->mK*elliptic_k(g->m2);
    double Tmo = g->mK*jacobi_icn(g->cos_i/sqrt(g->m2p),g->m2);
    double Pmp = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, 0.0);
    double Pmo = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, g->cos_i);
    double Pme = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, fabs(mu_e));
    // TODO it should be checked if that works for photons with m2m<0


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
 * /
    return time-RMAX;
    */
}




DEVICEFUNC
double geodesic_find_midplane_crossing(geodesic *g, int order)
// the value of T-integral in the equatorial plane
{
    if (g->q<=0.0) {
        #ifndef CUDA
        warning("WRN: q<0\n");
        #endif
        return NAN;
    }

    double u = g->cos_i/sqrt(g->m2p);
    if (!ensure_range(&u, -1.0, +1.0, 1e-4)) {
        #ifndef CUDA
        warning("u out of range (%e)\n", u);
        #endif
        return NAN;
    }

    // determine orientation of momentum vector (poloidal component of it)
    //if (dk2dx) (*dk2dx) = (order%2==0)?-1.:+1.;

    double pos;
    if (g->beta > 0.0)
        pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->m2) + jacobi_icn(u,g->m2) );
    else if (g->beta < 0.0)
        pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->m2) - jacobi_icn(u,g->m2) );
    else
        pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->m2) );

    if (pos > 2.*g->Rpa) pos = NAN;

    if (isnan(pos)) {
        #ifndef CUDA
        //warning("int_theta is nan (g.m2=%e)\n", g->m2);
        #endif
    }

    return pos;
}




DEVICEFUNC
void geodesic_follow(geodesic *g, double step, double *P, double *r, double *m, int *status)
// follows geodesic by evolution of g.x from its initial point (has to be set prior to call to follow())
// P, r, m has to be valid values on initial call
{
    const double MAXSTEP_FACTOR = 1e-2;
    double rbh = r_bh(g->a);

    do {
        double truestep = step/fabs(step) * min(fabs(step), MAXSTEP_FACTOR*sqrt(*r));
        (*P) = (*P) + truestep/(sqr(*r)+sqr((g->a)*(*m)));   // d(afp)/d(x) = r^2 + a^2*m^2
        (*r) = geodesic_position_rad(g, *P);
        (*m) = geodesic_position_pol(g, *P);
        if ((*r) < 1.01*rbh) {
            *status = 0;
            return;
        }
        step -= truestep;
    } while (fabs(step) > 1e-5);

    //int brp;
    //double ksgn_2;
    //photon_momentum(g->a, _r, _m, g->l, g->q, brp?-1:+1, ksgn_2, g->k);
    //void photon_momentum(double a, double r, double m, double l, double q2, double r_sign, double m_sign, double k[4])

    *status = 1;
}



































//------------------------------------------------------------------------------
// private methods
//------------------------------------------------------------------------------


DEVICEFUNC
int geodesic_priv_R_roots(geodesic *g, int *status)
// find roots of the R-integral
// http://adsabs.harvard.edu/abs/1998NewA....3..647C
{
    double a  = g->a;
    double l  = g->l;
    double q  = g->q;
    double a2 = sqr(a);
    double l2 = sqr(l);
    double A,B,C,D,E,F,X,Z,z;

    C = sqr(a-l)+q;
    D = 2./3.*(q+l2-a2);
    E = 9./4.*sqr(D)-12.*a2*q;
    F = -27./4.*sqr3(D)-108.*a2*q*D+108.*sqr(C);
    X = sqr(F)-4.*sqr3(E);
    if (X >= 0) {
        A = (F>sqrt(X)?+1:-1)*1./3.*pow(fabs(F-sqrt(X))/2.,1./3.)+ (F>-sqrt(X)?+1:-1)*1./3.*pow(fabs(F+sqrt(X))/2.,1./3.);
    } else {
         Z = sqrt(pow(F/54.,2) + pow(sqrt(-X)/54.,2));
         z = atan2(sqrt(-X)/54.,F/54.);
         A = pow(Z,1./3.)*2.*cos(z/3.);
    }
    B = sqrt(A+D);
    g->r1 = +B/2. + .5*csqrt(makeComplex(-A+2.*D-4.*C/B,0.0));
    g->r2 = +B/2. - .5*csqrt(makeComplex(-A+2.*D-4.*C/B,0.0));
    g->r3 = -B/2. + .5*csqrt(makeComplex(-A+2.*D+4.*C/B,0.0));
    g->r4 = -B/2. - .5*csqrt(makeComplex(-A+2.*D+4.*C/B,0.0));
    sort_roots(&g->nrr, &g->r1, &g->r2, &g->r3, &g->r4);

    if (g->nrr == 0) {
        // we do not yet know how to handle the case of non-equatorial photons
        // skip it silently
        if (status) *status = 0;
        #ifndef CUDA
        //fprintf(stderr,"nrr=0\n");
        #endif
        return FALSE;
    }

    // R-integral at turning point
    // value of position integral along geodesics at its periastron (radial turning point).
    // The function gives the value of the integral
    //     P = \int 1/\sqrt{R} dr = \int 1/\sqrt{\Theta} d\theta
    // This integral is integrated from infinity to its radial turning point.
    double r1,r2,r3,r4,u,v;
    if (g->nrr == 4) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        r3 = creal(g->r3);
        r4 = creal(g->r4);

        g->rp  = r1;
        if (r1==r2) {
            #ifndef CUDA
            fprintf(stderr, "ERR: no solution implemented for r1==r2\n");
            #endif
            g->Rpa = 0.0;
        } else {
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            g->Rpa = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r2-r4)/(r1-r4)), m4);
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
        g->rp  = r1;
        g->Rpa = 1./sqrt(A*B) * jacobi_icn((A-B)/(A+B), m2);
    } else {
        #ifndef CUDA
        fprintf(stderr, "ERR: no solution implemented for r1-r4 all complex\n");
        #endif
        g->rp  = 0.0;
        g->Rpa = 0.0;
    }

    return TRUE;
}




DEVICEFUNC
int geodesic_priv_T_roots(geodesic *g, int *status)
// T integral roots
// M(m) = q + (a2 - l2 - q)m2 - a2*m4 = -a2*(m4 + m2*(q + l2 -a2)/a2 - q/a2) = a2(m2m+m2)(m2p-m2)
{
    double a  = g->a;
    double l  = g->l;
    double q  = g->q;
    double a2 = sqr(a);
    double l2 = sqr(l);
    double qla, X;

    // the straightforward solution of quadratic function fails for small spins due to numerical cancelation:
    //   g->m2m = 1./(2.*a2) * ( sqrt(sqr(qla)+4.*q*a2) + qla );
    //   g->m2p = 1./(2.*a2) * ( sqrt(sqr(qla)+4.*q*a2) - qla );
    // so the roots mahe to be calculated with the help of equality m2m*m2p=q/a^2
    // note: Dexter&Agol define m2m with opposite sign
    qla = q + l2 - a2;
    X = sqrt(sqr(qla)+4.*q*a2) + qla;
    if (fabs(X)<1e-16) X=1e-16;
    g->m2m = X/2./a2;
    g->m2p = 2.*q/X;

    // for q>0.0 (photons that cross eq plane):
    // both m2p and m2m must be non-negative since m2p*m2m=q/a^2 > 0
    if ((q>0.0) && !ensure_range(&g->m2p, 0.0, 1.0, 1e-5)) {
        #ifndef CUDA
        fprintf(stderr,"m2p<0 (%e/%f/%f)\n", g->m2p,l,q);
        fprintf(stderr,"T err1\n");
        #endif
        if (status) *status = 0;
        return FALSE;
    }

    // for q<0.0 (photons that do not cross eq plane):
    if (q < 0.0) {
        // tbd
        if (status) *status = 0;
        return FALSE;
    }

    g->m2 = (g->m2m>0.0) ? g->m2p/(g->m2p+g->m2m) : g->m2p/(g->m2p-g->m2m);

    if (g->m2>=1.0) {
        *status=0;
        #ifndef CUDA
        fprintf(stderr,"WRN:g->m2>=1.0\n");
        fprintf(stderr,"T err2\n");
        #endif
        return FALSE;
    }
    if (!ensure_range(&g->m2, 0.0, 1.0, 1e-5)) {
        #ifndef CUDA
        fprintf(stderr,"m2<0 (%e)\n", g->m2);
        fprintf(stderr,"T err3\n");
        #endif
        if (status) *status = 0;
        return FALSE;
    }

    g->mK = 1./sqrt(a2*(g->m2p+g->m2m));

    //fprintf(stderr,"m2m=%.3e\n", g->m2m);
    //fprintf(stderr,"m2p=%.3e\n", g->m2p);
    //fprintf(stderr,"mK =%.3e\n", g->mK);

    return TRUE;
}









//------------------------------------------------------------------------------
// old routines (to be removed)
//------------------------------------------------------------------------------




DEVICEFUNC
double geodesic_s2i_int_R(geodesic *g, double r, int bpa)
//! Returns the value of position integral along geodesics at point r.
//!
//! The function gives the value of the integral
//!     P = \int 1/\sqrt{R} dr = \int 1/\sqrt{\Theta} d\theta
//! This integral is integrated from infinity to the point, where the geodesic
//! reaches radius r either before or behind its periastron.
//! The value of the integral increases monotonicly from infinity along the geodesic
//! and thus we use it for parametrizing it.
//! Note it is not affine parameter, which would be other choice for parametrization.
//!
//! Parameters:
//!     g      - geodesic
//!     r      - radial coordinate [GM/c^2]
//!     bpa    - position with respect to periastron (0=before periastron, 1=after periastron)
//!
//! Returns:
//!     Value of the position integral between infinity and given point.
{
    double r1,r2,r3,r4,u,v;

    #ifdef CUDA
    if (r  < g->rp) asm("exit;");
    #else
    if (r  < g->rp) error("geodesic_s2i_int_R: r < periastron (%.3e < %.3e)", r, g->rp);
    #endif
    if (r == g->rp) return g->Rpa;

    if (g->nrr == 4) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        r3 = creal(g->r3);
        r4 = creal(g->r4);

        if (r1==r2) {
            #ifndef CUDA
            fprintf(stderr, "ERR (geodesic_affine): no solution implemented for r1==r2\n");
            #endif
            return 0.0;
        } else {
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            double R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt(((r2-r4)*(r-r1))/((r1-r4)*(r-r2))), m4);
            return (bpa) ? g->Rpa + R : g->Rpa - R;
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
        return (bpa) ? g->Rpa + R : g->Rpa - R;
    } else {
        #ifndef CUDA
        fprintf(stderr, "ERR: no solution implemented for r1-r4 all complex\n");
        #endif
        return 0.0;
    }
}

DEVICEFUNC
double geodesic_s2i_inv_R(geodesic *g, double r, int* bpa)
//! Returns the value of position integral along geodesics at point r.
//!
//! The function gives the value of the integral
//!     P = \int 1/\sqrt{R} dr = \int 1/\sqrt{\Theta} d\theta
//! This integral is integrated from infinity to the point, where the geodesic
//! reaches radius r either before or behind its periastron.
//! The value of the integral increases monotonicly from infinity along the geodesic
//! and thus we use it for parametrizing it.
//! Note it is not affine parameter, which would be other choice for parametrization.
//!
//! Parameters:
//!     g      - geodesic
//!     r      - radial coordinate [GM/c^2]
//!     bpa    - position with respect to periastron (0=before periastron, 1=after periastron)
//!
//! Returns:
//!     Value of the position integral between infinity and given point.
{
    double r1,r2,r3,r4,u,v;

    #ifdef CUDA
    if (r  < g->rp) asm("exit;");
    #else
    if (r  < g->rp) error("geodesic_s2i_int_R: r < periastron (%.3e < %.3e)", r, g->rp);
    #endif
    if (r == g->rp) return g->Rpa;

    if (g->nrr == 4) {
        r1 = creal(g->r1);
        r2 = creal(g->r2);
        r3 = creal(g->r3);
        r4 = creal(g->r4);

        if (r1==r2) {
            #ifndef CUDA
            fprintf(stderr, "ERR (geodesic_affine): no solution implemented for r1==r2\n");
            #endif
            return 0.0;
        } else {
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            double R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt(((r2-r4)*(r-r1))/((r1-r4)*(r-r2))), m4);
            return (bpa) ? g->Rpa + R : g->Rpa - R;
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
        return (bpa) ? g->Rpa + R : g->Rpa - R;
    } else {
        #ifndef CUDA
        fprintf(stderr, "ERR: no solution implemented for r1-r4 all complex\n");
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
    if (!ensure_range(&u, -1.0, +1.0, 1e-4)) {
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
    //fprintf(stderr,"intT: %d K=%.2e psi=%.2e u=%.2e res=%.2e rpr=%.2e/%.2e\n", order, g->mK*elliptic_k(g->m2), sgn_psi*g->mK*psi, u, res, g->Rpa, g->rp);

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
    #ifndef CUDA
    error("geodesic_s2i_inv_T: ???");
    #endif
    return 0.0;
}


DEVICEFUNC
double geodesic_s2i_timedelay(geodesic *g, double x, double *opt_r, double *opt_m)
{
    int    bpr  = ((g->nrr==4) && (x > g->Rpa));              // beyond periastron radius
    double a2 = sqr(g->a);
    double rp = 1. + sqrt(1.-a2);
    double rm = 1. - sqrt(1.-a2);
    double T  = 0.0;
    double r    = (opt_r) ? *opt_r : geodesic_s2i_inv_R(g,x,&bpr);
    //unused double mu_e = (opt_m) ? *opt_m : geodesic_s2i_inv_T(g,x,NULL);

    double RMAX = 1000.;
    #ifdef CUDA
    if (r > RMAX) asm("exit;");
    if (r < g->rp) asm("exit;");
    #else
    if (r > RMAX) error("geodesic_s2i_timedelay: r > RMAX");
    if (r < g->rp) error("geodesic_s2i_timedelay: r < r_p");
    #endif

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
    int    bpr  = (g->nrr==4) && (x > g->Rpa);              // beyond periastron radius
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
        #ifdef CUDA
        if (bpr) asm("exit;");
        #else
        if (bpr) error("geodesic_s2i_phi: cannot be (beyond_pa) && (nrr==2)");
        #endif
        double r1 = creal(g->r1);
        double r2 = creal(g->r2);
        double A = integral_R_rp_cc2_inf(r1, r2, g->r3, rp, r);
        double B = integral_R_rp_cc2_inf(r1, r2, g->r3, rm, r);
        phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
    } else {
        #ifdef CUDA
        asm("exit;");
        #else
        error("geodesic_s2i_phi: g->nrr != [2,4]");
        #endif
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
    int    bpr  = (g->nrr==4) && (x > g->Rpa);              // beyond periastron radius
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
        #ifndef CUDA
        error("geodesic_s2i_afp: g->nrr != [2,4]");
        #endif
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
DEVICEFUNC
void geodesic_s2i_solution_surface(geodesic *g, double *r, double *m, double (*H)(double), double accuracy, int *status)
{
    *r = 0.0; *m = 0.0;
    double rh = r_bh(g->a);
    double r0 = max(1000.0, g->rp*2.);
    double r1, m1, R1, h1, H1;
    double r2, m2;//, x2, R2, h2, H2;

    do {
      // initialize photon-follow from r0
      geodesic_s2i_follow_init(g, r0, status);
      if (!(*status)) return;

      // get initial position
      //x1 = g->x;
      r1 = r2 = geodesic_s2i_inv_R(g, g->x, NULL);
      m1 = m2 = geodesic_s2i_inv_T(g, g->x, NULL);
      R1 = r1*sqrt(1.-m1*m1);
      h1 = r1*m1;
      H1 = fabs((*H)(R1));

      r0 *= 2.0;
    } while ((h1<=H1) && (r0<5e4));

    // test if the initial point is already under surface
    if (h1 <= H1) { *status=0; return; }

    // follow trajectory with adaptive step until surface is reached
    do {
        double step;
        r2=r1; m2=m1; //R2=R1; x2=x1; h2=h1; H2=H1;

        // new step
        step=max((h1-H1)/5., H1*accuracy);
        geodesic_s2i_follow(g, step, &r1, &m1, NULL, NULL, status);
        //fprintf(stderr,"r=%.3f  step=%.3e H1=%.3e, sta=%d\n", r1, step, H1, *status);
    	if (!*status) return;

        // current photon and surface location
        //x1 = g->x;
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
void x_geodesic_s2i_follow_init(geodesic *g, double rmax, int *status)
{
/ *
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
* /
    #ifndef CUDA
    if (rmax<10.0) fprintf(stderr, "WRN: rmax should be >10 (geodesic_s2i_follow_init)\n");
    #endif
    *status = 0;
    if (rmax < g->rp) return;
    g->x = geodesic_s2i_int_R(g, rmax, 0);
    *status = 1;
}


DEVICEFUNC
void x_geodesic_s2i_follow(geodesic *g, double step, double *r, double *m, double *f, double *t, int *status)
{
/ *
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
* /

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
*/

#undef theta_int
#undef theta_inv
