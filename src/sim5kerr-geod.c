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
#define theta_int(x) (g->mK*jacobi_icn((x)/sqrt(g->m2p),g->mm))
#define theta_inv(x) (sqrt(g->m2p)*jacobi_cn((x)/g->mK,g->mm))

// unit-private function declarations
DEVICEFUNC double geodesic_priv_RR(geodesic *g, double r);
DEVICEFUNC double geodesic_priv_TT(geodesic *g, double m);
DEVICEFUNC int geodesic_priv_R_roots(geodesic *g, double r0, int *error);
DEVICEFUNC int geodesic_priv_T_roots(geodesic *g, double m0, int *error);



DEVICEFUNC
int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *error)
//! Makes setup for geodesic that is specified by impact parameters at infinity.
//!
//! Parameters:
//!     i      - inclination angle of observer (angle between BH rotation axis and
//!              direction to observer) [radians]
//!     a      - BH spin [0..1]
//!     alpha  - impact parameter in horizontal direction [GM/c^2]
//!     beta   - impact parameter in vertical direction [GM/c^2]
//!     g      - structure with information about geodesic
//!     error  - error code
//!
//! Returns:
//!     * error is 0 if setup is sucessfull or contains a non-zero error code
//!     * information about geodesic is stored in structure g
{

    if ((a < 0.0) || (a > 1.-1e-6)) {
        if (error) *error = GD_ERROR_SPIN_RANGE;
        return FALSE;
    }

    if ((i <= 0.0) || (i>=PI_half)) {
        if (error) *error = GD_ERROR_INCL_RANGE;
        return FALSE;
    }

    if (beta == 0.0) beta = +1e-6;

    g->a = max(1e-8, a);
    g->incl  = i;
    g->cos_i = cos(i);
    g->alpha = alpha;
    g->beta  = beta;

    // constants of motion
    g->l = -alpha*sin(i);
    g->q = sqr(beta) + sqr(cos(i))*(sqr(alpha)-sqr(a));

    if (g->q == 0.0) {
        // q=0 are trajectories confined to the equatorial plane
        // we are not yet ready to handle this case (Dexter&Agol provide a solution in case of need)
        if (error) *error = GD_ERROR_Q_RANGE;
        return FALSE;
    }

    // get geodesic
    if (!geodesic_priv_R_roots(g, DBL_MAX, error)) return FALSE;
    if (!geodesic_priv_T_roots(g, g->cos_i, error)) return FALSE;

    // value of T-integral between turning points \int[-\mu_plus..\mu_plus]
	g->Tpp = 2.*theta_int(0.0);

    // value of T-integral between observer's position and turning point \int[cos_i..\mu_plus]
	g->Tip = theta_int(g->cos_i);

    if (error) *error = GD_OK;
    return TRUE;
}




DEVICEFUNC
int geodesic_init_src(double a, double r, double m, double k[4], int ppc, geodesic *g, int *error)
//! Makes setup for geodesic that is specified by a point and direction (4-momentum vector).
//!
//! Parameters:
//!     a      - BH spin [0..1]
//!     r      - radial coordinate [GM/c^2]
//!     m      - poloidal coordinate [cos(theta)]
//!     k      - direction of the photon (4-momentum vector)
//!     ppc    - position with respect to periastron (0=before periastron, 1=after periastron)
//!     g      - structure with information about geodesic
//!     error  - error code [out]
//!
//! Returns:
//!     * error is GD_OK if setup is sucessfull, in case of an error it contains a non-zero error code
//!     * information about geodesic is stored in structure g
{
    // calculate motion constants
    double l,q;
    photon_motion_constants(a, r, m, k, &l, &q);

    g->a = max(1e-8, a);
    g->l = l;
    g->q = q;

    // ...
    g->cos_i = g->alpha = g->beta = NAN;

    if (!geodesic_priv_R_roots(g, r, error)) return FALSE;
    if (!geodesic_priv_T_roots(g, m, error)) return FALSE;
    
    // r is bellow radial turning point - this geodesic cannot escape
    //if (r < g->rp) {
    //    //TODO: should be implemented at some point
    //    if (error) *error = GD_ERROR_BOUND_GEODESIC;
    //    return FALSE;
    //}

    // determine at-infinity parameters of the geodesic (inclination, impact parameters)
    // - inclination is poloidal coordinate to which the photon in its current direction will arrive
    if (isnan(g->cos_i) && (r > g->rp)) {
        double T, Tmp, Tpp, sign_dm;
        Tmp = theta_int(m);                         // int_{\mu}^{\mu_plus}
        Tpp = 2.*theta_int(0.0);                    // int_{-\mu_plus}^{\mu_plus}
        T = geodesic_P_int(g, r, ppc);              // int_{r}^{\inf}
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
        g->incl  = acos(g->cos_i);
        g->alpha = -g->l/sqrt(1.0-sqr(g->cos_i));
        g->beta  = -sign_dm * sqrt(g->q - sqr(g->cos_i)*(sqr(g->alpha)-sqr(g->a)));
    }

    // value of T-integral between turning points \int[-\mu_plus..+\mu_plus]
	g->Tpp = 2.*theta_int(0.0);

    // value of T-integral between observer's position and turning point \int[cos_i..\mu_plus]
	g->Tip = theta_int(g->cos_i);

    if (error) *error = GD_OK;
    return TRUE;
}




DEVICEFUNC
double geodesic_P_int(geodesic *g, double r, int ppc)
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
//!     ppc    - position with respect to periastron (0=before periastron, 1=after periastron)
//!
//! Returns:
//!     Value of the position integral between infinity and given point.
{
    double r1,r2,r3,r4,u,v, mm, R, A, B;

    #ifdef CUDA
    if (r  < g->rp) asm("exit;");
    #else
    if (r  < g->rp) error("(geodesic_P_int): r < periastron (%.3e < %.3e; nrr=%d)", r, g->rp, g->nrr);
    #endif
    if (r == g->rp) return g->Rpc;


    switch (g->class) {
        case CLASS_RR:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt(((r2-r4)*(r-r1))/((r1-r4)*(r-r2))), mm);
            return (ppc) ? g->Rpc + R : g->Rpc - R;
        
        case CLASS_RR_DBL:
            #ifndef CUDA
            error("(geodesic_P_int): not implemented for CLASS_RR_DBL");
            #endif
            return NAN;

        case CLASS_RR_BH:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r1-r3)/(r2-r3)*(r2-r)/(r1-r)), mm);
            return (ppc) ? g->Rpc + R : g->Rpc - R;

        case CLASS_RC:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            u  = creal(g->r3);
            v  = cimag(g->r3);
            A = sqrt(sqr(r1-u)+sqr(v));
            B = sqrt(sqr(r2-u)+sqr(v));
            mm = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
            R  = 1./sqrt(A*B) * jacobi_icn(((A-B)*r+r1*B-r2*A)/((A+B)*r-r1*B-r2*A), mm);
            //if (R > g->Rpc) fprintf(stderr, "RC: res=%e Rpc=%e R=%e r=%.4e r1=%.4e r2=%.4e A=%.4e B=%.4e mm=%.4e z=%.4e icn=%.4e\n", g->Rpc-R, g->Rpc, R, r, r1, r2, A, B, mm, ((A-B)*r+r1*B-r2*A)/((A+B)*r-r1*B-r2*A), jacobi_icn(((A-B)*r+r1*B-r2*A)/((A+B)*r-r1*B-r2*A), mm));
            // RC case has no turning point, so P can be 0 < P < Rpc
            return g->Rpc-R;

        case CLASS_CC:
            r1 = creal(g->r1);  //b1
            r2 = creal(g->r3);  //b2
            r3 = cimag(g->r1);  //a1
            r4 = cimag(g->r3);  //a2
            A = sqrt(sqr(r1-r2) + sqr(r3+r4));
            B = sqrt(sqr(r1-r2) + sqr(r3-r4));
            double g1 = sqrt((4.*sqr(r3)-sqr(A-B))/(sqr(A+B)-4.*sqr(r3)));
            mm = 4.*A*B/sqr(A+B);
            R = 2./(A+B) * jacobi_itn((r-r1+r3*g1)/(r3+r1*g1-g1*r), mm);
            if ((R)<0) {
                fprintf(stderr, "CC: res=%e ppc=%d Rpc=%e R=%e\n", R, ppc, g->Rpc, R);
            }
            return g->Rpc-R;
    }

    // this point should never be reached
    return NAN;
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

    if ((P<=0.0)||(P>=2.*g->Rpc)) {
        #ifndef CUDA
        error("(geodesic_position_rad) P out of range (P=%e, 2Rpa=%e)", P, 2*g->Rpc);
        #endif
        return NAN;
    }
    if (P == g->Rpc) return g->rp;

    switch (g->class) {
        case CLASS_RR:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            double x4 = 0.5*fabs    (P - g->Rpc)*sqrt((r2-r4)*(r1-r3));
            double sn2 = pow( jacobi_sn(x4,m4), 2.0);
            return ( r1*(r2-r4)-r2*(r1-r4)*sn2 ) / ( r2-r4-(r1-r4)*sn2 );

        case CLASS_RR_DBL:
            #ifndef CUDA
            error("(geodesic_position_rad): not implemented for CLASS_RR_DBL");
            #endif
            return NAN;

        case CLASS_RR_BH:
            #ifndef CUDA
            error("(geodesic_position_rad): not implemented for CLASS_RR_BH");
            #endif
            return NAN;

        case CLASS_RC:
            // RC case has no turning point, so P cannot be larger than Rpc
            if (P > g->Rpc) return NAN;
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            u  = creal(g->r3);
            v  = cimag(g->r3);
            double A = sqrt(sqr(r1-u)+sqr(v));
            double B = sqrt(sqr(r2-u)+sqr(v));
            double m2 = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
            double cn = jacobi_cn(sqrt(A*B)*(g->Rpc-P), m2);
            fprintf(stderr, "r-RC: mm=%.4e z=%.4e icn=%.4e  r=%.4e\n", m2, sqrt(A*B)*(g->Rpc-P), cn, (r2*A - r1*B - (r2*A+r1*B)*cn ) / ( (A-B) - (A+B)*cn ));
            return (r2*A - r1*B - (r2*A+r1*B)*cn ) / ( (A-B) - (A+B)*cn );
        
        case CLASS_CC:
            #ifndef CUDA
            error("(geodesic_position_rad): not implemented for CLASS_CC");
            #endif
            return NAN;
    }

    // this point should never be reached
    return NAN;
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
    double sign_dm, T;

    switch (g->class) {
        // trajectories that go to infinity
        case CLASS_RR:
        case CLASS_RC:
        case CLASS_CC:
            // sign_dm = d(m)/d(P)
            // at infinity, d(m)/d(P) > 0 if beta>0, and d(m)/d(P) < 0 if beta<0
            sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
            T = (sign_dm>0.0) ? -(g->Tpp-g->Tip) : -(g->Tip);
            // d(m)/d(P) changes sign each multiple of Tpp=\int_{-mu_+}^{+mu_+}
            while (P > T+g->Tpp) {
                T += g->Tpp;
                sign_dm = -sign_dm;
                //fprintf(stderr,"P=%.3e T=%.3e T+Tpp=%.3e delta=%.3e\n", P, T, T+g->Tpp, P-T-g->Tpp);
            }
            return -sign_dm*theta_inv(P-T);

        case CLASS_RR_DBL:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for CLASS_RR_DBL");
            #endif
            return NAN;

        case CLASS_RR_BH:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for CLASS_BH");
            #endif
            return NAN;
    }

    // default
    return NAN;
}




DEVICEFUNC
double geodesic_position_pol_sign_k_theta(geodesic *g, double P)
//! Gives the sign of k^\theta component of the 4-momentum.
//!
//! Parameters:
//!     g      - geodesic
//!     P      - value of the position integral
//!
//! Returns:
//!     +1 or -1
{
    double sign_dm, T;

    switch (g->class) {
        // trajectories that go to infinity
        case CLASS_RR:
        case CLASS_RC:
        case CLASS_CC:
            // sign_dm = d(m)/d(P)
            // at infinity, d(m)/d(P) > 0 if beta>0, and d(m)/d(P) < 0 if beta<0
            sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
            T = (sign_dm>0.0) ? -(g->Tpp-g->Tip) : -(g->Tip);
            // d(m)/d(P) changes sign each multiple of Tpp=\int_{-mu_+}^{+mu_+}
            while (P > T+g->Tpp) {
                T += g->Tpp;
                sign_dm = -sign_dm;
            }
            // dk[2] = -d(m)
            return (sign_dm<0) ? +1 : -1;

        case CLASS_RR_DBL:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for CLASS_RR_DBL");
            #endif
            return NAN;

        case CLASS_RR_BH:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for CLASS_BH");
            #endif
            return NAN;
    }

    // default
    return NAN;
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

    int    ppc  = (g->nrr>0) && (P > g->Rpc);              // beyond periastron
    double a2 = sqr(g->a);
    double rp   = 1. + sqrt(1.-a2);
    double rm   = 1. - sqrt(1.-a2);
    double r1, r2, r3, r4, A, B;

    switch (g->class) {
        // trajectories that go to infinity
        case CLASS_RR:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            A = integral_R_rp_re_inf(r1, r2, r3, r4, rp) + (ppc?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rp, r);
            B = integral_R_rp_re_inf(r1, r2, r3, r4, rm) + (ppc?+1:-1)*integral_R_rp_re(r1, r2, r3, r4, rm, r);
            phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
            break;

        case CLASS_RR_DBL:
            #ifndef CUDA
            error("(geodesic_position_azm): not implemented for CLASS_RR_DBL");
            #endif
            return NAN;

        case CLASS_RR_BH:
            #ifndef CUDA
            error("(geodesic_position_azm): not implemented for CLASS_RR_BH");
            #endif
            return NAN;

        case CLASS_RC:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            A = integral_R_rp_cc2_inf(r1, r2, g->r3, rp, r);
            B = integral_R_rp_cc2_inf(r1, r2, g->r3, rm, r);
            phi += 1./sqrt(1.-a2) * ( A*(g->a*rp-g->l*a2/2.) - B*(g->a*rm-g->l*a2/2.) );
            //if (isnan(phi)) fprintf(stderr,"phi nan already 2 (%e %e)\n",g->a,rm);
            break;

        case CLASS_CC:
            #ifndef CUDA
            error("(geodesic_position_azm): not implemented for CLASS_CC");
            #endif
            return NAN;
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
double geodesic_timedelay(geodesic *g, double P1, double r1, double m1, double P2, double r2, double m2)
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
    double time = 0.0;

    // exchange endpoints if needed, so it can be assumed P1 < P2
    if (P1 > P2) {
        double tmp;
        tmp=P2; P2=P1; P1=tmp;
        tmp=r2; r2=r1; r1=tmp;
        tmp=m2; m2=m1; m1=tmp;
    }

    if (r1 == 0) {
        r1 = geodesic_position_rad(g, P1);
        m1 = geodesic_position_pol(g, P1);
    }

    if (r2 == 0) {
        r2 = geodesic_position_rad(g, P2);
        m2 = geodesic_position_pol(g, P2);
    }


    double a2 = sqr(g->a);
    double rp   = 1. + sqrt(1.-a2);
    double rm   = 1. - sqrt(1.-a2);
    double ra = creal(g->r1);
    double rb = creal(g->r2);
    double rc = creal(g->r3);
    double rd = creal(g->r4);
    double R0, R1, R2, RA, RB, A, B, s;

    //fprintf(stderr,"l=%e q=%e\n", g->l, g->q);
    //fprintf(stderr,"r1=%e m1=%e P1=%e\n", r1, m1, P1);
    //fprintf(stderr,"r2=%e m2=%e P2=%e\n", r2, m2, P2);

    if (r1 < g->rp) error("geodesic_timedelay: r1 < r_p (%e/%e)", r1, g->rp);
    if (r2 < g->rp) error("geodesic_timedelay: r2 < r_p (%e/%e)", r2, g->rp);

    switch (g->class) {
        // trajectories that go to infinity
        case CLASS_RR:
            s = (((P1 > g->Rpc)&&(P2 < g->Rpc)) || ((P1 < g->Rpc)&&(P2 > g->Rpc))) ? +1 : -1;
            R0 = integral_R_r0_re(ra, rb, rc, rd, r1)     + s*integral_R_r0_re(ra, rb, rc, rd, r2);
            R1 = integral_R_r1_re(ra, rb, rc, rd, r1)     + s*integral_R_r1_re(ra, rb, rc, rd, r2);
            R2 = integral_R_r2_re(ra, rb, rc, rd, r1)     + s*integral_R_r2_re(ra, rb, rc, rd, r2);
            RA = integral_R_rp_re(ra, rb, rc, rd, rp, r1) + s*integral_R_rp_re(ra, rb, rc, rd, rp, r2);
            RB = integral_R_rp_re(ra, rb, rc, rd, rm, r1) + s*integral_R_rp_re(ra, rb, rc, rd, rm, r2);
            A = (-g->a*g->l+4.)*rp - 2.*a2;
            B = (+g->a*g->l-4.)*rm + 2.*a2;
            time += 4.*fabs(R0) + 2.*fabs(R1) + fabs(R2) + (A*fabs(RA) + B*fabs(RB))/sqrt(1.-a2);
            //fprintf(stderr,"RR1=%e  RR2=%e\n", geodesic_priv_RR(g,r1), geodesic_priv_RR(g,r2));
            //fprintf(stderr,"ra=%.3e rb=%.3e rc=%.3e rd=%.3e\n", ra, rb, rc, rd);
            //fprintf(stderr,"R0=%.3e R1=%.3e R2=%.3e RA=%.3e RB=%.3e A=%.3e B=%.3e dt=%.3e \n", R0,R1,R2,RA,RB,A,B,time);
            break;

        case CLASS_RR_DBL:
            #ifndef CUDA
            error("(geodesic_timedelay): not implemented for CLASS_RR_DBL");
            #endif
            return NAN;

        case CLASS_RR_BH:
            #ifndef CUDA
            error("(geodesic_timedelay): not implemented for CLASS_RR_BH");
            #endif
            return NAN;

        case CLASS_RC:
            R0 = integral_R_r0_cc(ra, rb, g->r3, r1) -integral_R_r0_cc(ra, rb, g->r3, r2);
            R1 = (r1<r2) ? integral_R_r1_cc(ra, rb, g->r3, r1, r2) : integral_R_r1_cc(ra, rb, g->r3, r2, r1);
            R2 = (r1<r2) ? integral_R_r2_cc(ra, rb, g->r3, r1, r2) : integral_R_r2_cc(ra, rb, g->r3, r2, r1);
            RA = (r1<r2) ? integral_R_rp_cc2(ra, rb, g->r3, rp, r1, r2) : integral_R_rp_cc2(ra, rb, g->r3, rp, r2, r1);
            RB = (r1<r2) ? integral_R_rp_cc2(ra, rb, g->r3, rm, r1, r2) : integral_R_rp_cc2(ra, rb, g->r3, rm, r2, r1);
            A = (-g->a*g->l+4.)*rp - 2.*a2;
            B = (+g->a*g->l-4.)*rm + 2.*a2;
            time += 4.*fabs(R0) + 2.*fabs(R1) + fabs(R2) + (A*fabs(RA) + B*fabs(RB))/sqrt(1.-a2);
            //fprintf(stderr,"RR1=%e  RR2=%e\n", geodesic_priv_RR(g,r1), geodesic_priv_RR(g,r2));
            //fprintf(stderr,"ra=%.3e rb=%.3e rc=%.3e rd=%.3e\n", ra, rb, rc, rd);
            //fprintf(stderr,"R0=%.3e R1=%.3e R2=%.3e RA=%.3e RB=%.3e A=%.3e B=%.3e dt=%.3e \n", R0,R1,R2,RA,RB,A,B,time);
            break;

        case CLASS_CC:
            #ifndef CUDA
            error("(geodesic_timedelay): not implemented for CLASS_CC");
            #endif
            return NAN;
    }

    /*

    // part 2 - T integral
    // TODO it should be checked if that works for photons with m2m<0
    double time_pp = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, 0.0)*2;        // M-integral from -mu_plus to +mu_plus
    double time_m1 = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, fabs(m1));     // M-integral from m1 to +mu_plus
    double time_m2 = sqr(g->a)*integral_T_m2(g->m2m, g->m2p, fabs(m2));     // M-integral from m2 to +mu_plus

    double dP = P2-P1;

    while (T >= g->Tpp) {
        T    -= g->Tpp;
        time += time_pp;
    }
    
    if (T > g->Tpp/2) {
        time += fabs(time_m2 - time_m1);
    } else {
        time += fabs(time_m2 - time_m1);
    }
    
    
    double sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
    if (sign_dm > 0.0) {
        T = -(g->Tpp-g->Tip);
        time -= time_pp-time_ip;
    } else {
        T = -g->Tip;
        time -= time_ip;
    }

    while (fabs(P2-P1) >= T+g->Tpp) {
        T += g->Tpp;
        time += time_pp;
        sign_dm = -sign_dm;
        break;
    }

    time += (sign_dm<0) ? time_mp : time_pp-time_mp;
*/
    return time;

/*
old:
    double Tmp = g->mK*elliptic_k(g->mm);
    double Tmo = g->mK*jacobi_icn(g->cos_i/sqrt(g->m2p),g->mm);
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
 */
    return time;
}




DEVICEFUNC
double geodesic_dm_sign(geodesic *g, double P)
//! Gives the sign of the derivative d(m)/d(P) at current position
//!
//! Parameters:
//!     g      - geodesic
//!     P      - value of the position integral
//!
//! Returns:
//!     Sign of d(m)/d(P), i.e. +1 or -1.
{
    double sign_dm, T;

    switch (g->class) {
        // trajectories that go to infinity
        case CLASS_RR:
        case CLASS_RC:
        case CLASS_CC:
            // sign_dm = d(m)/d(P)
            // at infinity, d(m)/d(P) > 0 if beta>0, and d(m)/d(P) < 0 if beta<0
            sign_dm = (g->beta>=0.0) ? +1.0 : -1.0;
            T = (sign_dm>0.0) ? -(g->Tpp-g->Tip) : -(g->Tip);

            // d(m)/d(P) changes sign each multiple of Tpp=\int_{-mu_+}^{+mu_+}
            while (P > T+g->Tpp) {
                T += g->Tpp;
                sign_dm = -sign_dm;
            }
            return sign_dm;

        case CLASS_RR_DBL:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for CLASS_RR_DBL");
            #endif
            return NAN;

        case CLASS_RR_BH:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for CLASS_BH");
            #endif
            return NAN;
    }

    // default
    return NAN;
}




DEVICEFUNC
void geodesic_momentum(geodesic *g, double P, double r, double m, double k[])
//! Gives the 4-momentum of photons at given position along the geodesic.
//! The function needs to know [r,m] coordinates of the point at the trajectory.
//! If both r=m=0.0, the required values are computed from the value P of 
//! the position integral. To save these computations, value of [r,m] coordinates
//! can be given to the function, if they have been computed before.
//! Note: It is important to give the correct values of [r,m] corresponding to current position.
//! The orientation of the momentum vector is always in the direction of increasing P,
//! i.e. it points towards the radial turning point before it is reached and away 
//! from the radial turning point after it is reached.
//!
//! Parameters:
//!     g      - geodesic
//!     P      - value of the position integral
//!     r      - radial coordinate (value or zero)
//!     m      - poloidal coordinate (value or zero)
//!     k      - 4-momentum vector (output)
//!
//! Returns:
//!     Photon 4-momentum.
{
    double dm;

    // calc [r,m] coordinate if it has not been provided
    if ((r==0.0) && (m==0.0)) {
        r = geodesic_position_rad(g, P);
        m = geodesic_position_pol(g, P);
    }

    switch (g->class) {
        // trajectories that go to infinity
        case CLASS_RR:
        case CLASS_RC:
        case CLASS_CC:
            dm = geodesic_dm_sign(g, P);
            photon_momentum(g->a, r, m, g->l, g->q, (P<g->Rpc?-1:+1), dm, k);
            return;
            
        case CLASS_RR_DBL:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for CLASS_RR_DBL");
            #endif
            k[0]=k[1]=k[2]=k[3]=NAN;
            return;

        case CLASS_RR_BH:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for CLASS_BH");
            #endif
            k[0]=k[1]=k[2]=k[3]=NAN;
            return;
    }

    return;
}




DEVICEFUNC
double geodesic_find_midplane_crossing(geodesic *g, int order)
// the value of T-integral in the equatorial plane
{
    if (g->q<=0.0) {
        // there is no midplane crossing for photons with q<=0
        return NAN;
    }

    double u = g->cos_i/sqrt(g->m2p);
    if (!ensure_range(&u, -1.0, +1.0, 1e-4)) {
        #ifndef CUDA
        error("(geodesic_find_midplane_crossing): u out of range (%e)", u);
        #endif
        return NAN;
    }

    double pos;
    if (g->beta > 0.0)
        pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->mm) + jacobi_icn(u,g->mm) );
    else if (g->beta < 0.0)
        pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->mm) - jacobi_icn(u,g->mm) );
    else
        pos = g->mK*( (2.*(double)order+1.)*elliptic_k(g->mm) );

    if (pos > 2.*g->Rpc) pos = NAN;

    return pos;
}




DEVICEFUNC
void geodesic_follow(geodesic *g, double step, double *P, double *r, double *m, int *status)
// follows geodesic by evolution of g.x from its initial point (has to be set prior to call to follow())
// P, r, m has to be valid values on initial call
{
    const double MAXSTEP_FACTOR = 5e-2;

    do {
        double truestep = step/fabs(step) * min(fabs(step), MAXSTEP_FACTOR*sqrt(*r));
        (*P) = (*P) + truestep/(sqr(*r)+sqr((g->a)*(*m)));   // d(afp)/d(x) = r^2 + a^2*m^2
        (*r) = geodesic_position_rad(g, *P);
        (*m) = geodesic_position_pol(g, *P);
        if ((*r) < 1.01*r_bh(g->a)) {
            if (status) *status = 0;
            return;
        }
        if ((*P < 0.0) || (*P>2.*g->Rpc)) {
            if (status) *status = 0;
            return;
        }
        step -= truestep;
    } while (fabs(step) > 1e-5);

    if (status) *status = 1;
}



































//------------------------------------------------------------------------------
// private methods
//------------------------------------------------------------------------------


DEVICEFUNC
double geodesic_priv_RR(geodesic *g, double r) 
//! Value of R(r) function of the R-integral at given R
//! Cadez,Fanton,Calvani(1998) Eq.(7)
{
    return sqr4(r) + (sqr(g->a) - sqr(g->l) - g->q)*sqr(r) + 2.*(g->q + sqr(g->l-g->a))*r - sqr(g->a)*g->q;
}


DEVICEFUNC
double geodesic_priv_TT(geodesic *g, double m) 
//! Value of Theta(m) function of the theta-integral at given m=cos(theta)
//! Cadez,Fanton,Calvani(1998) Eq.(10)
{
    return sqr(g->a)*(g->m2m + m*m)*(g->m2p - m*m);
}


DEVICEFUNC
int geodesic_priv_R_roots(geodesic *g, double r0, int *error)
// find roots of the R-integral
// http://adsabs.harvard.edu/abs/1998NewA....3..647C (http://web.oapd.inaf.it/calvani/cadez.ps)
{
    double a  = g->a;
    double l  = g->l;
    double q  = g->q;
    double a2 = sqr(a);
    double l2 = sqr(l);
    double A,B,C,D,E,F,X,Z,z;

    // calculate roots of R(r) after Cadez et al (1998)
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

    // trajectory class
    switch (g->nrr) {
        case 4:
            g->class = CLASS_RR;
            // r0 can only be between r3 and r2; or it can be above r1
            // anything else is an error
            if ((r0<creal(g->r3)) || ((r0>creal(g->r2)) && (r0<creal(g->r1)))) {
                if (error) *error = GD_ERROR_UNKNOWN_SOLUTION;
                return FALSE;
            }
            // if r1 and r2 are close to each other, it is a double root solution
            if (fabs(creal(g->r1)-creal(g->r2)) < 1e-8) {
                g->class = CLASS_RR_DBL;
                if (error) *error = GD_ERROR_CLASS_RR_DOUBLE;
                return FALSE;
            }
            // if r0 is between r3 and r2, it is an inner solution
            if ((r0>=creal(g->r3)) && (r0<=creal(g->r2))) {
                g->class = CLASS_RR_BH;
            }
            break;
        case 2:
            g->class = CLASS_RC;
            break;
        case 0:
            g->class = CLASS_CC;
            break;
        default:
            if (error) *error = GD_ERROR_UNKNOWN_SOLUTION;
            return FALSE;
    }

    // set radius of turning point and integral value at turning point
    double r1,r2,r3,r4,u,v,mm;
    switch (g->class) {
        case CLASS_RR:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            g->rp  = r1;
            g->Rpc = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r2-r4)/(r1-r4)), mm);
            break;

        case CLASS_RR_BH:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            g->rp  = r2;
            g->Rpc = 2./sqrt((r1-r3)*(r2-r4)) * elliptic_k(mm); //ellipticK(mm) = jacobi_isn(1, mm)
            break;

        case CLASS_RC:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            u  = creal(g->r3);
            v  = cimag(g->r3);
            A = sqrt(sqr(r1-u)+sqr(v));
            B = sqrt(sqr(r2-u)+sqr(v));
            mm = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
            g->rp  = r1;
            g->Rpc = 1./sqrt(A*B) * jacobi_icn((A-B)/(A+B), mm);
            break;

        case CLASS_CC:
            r1 = creal(g->r1);   //b1
            r2 = creal(g->r3);   //b2
            r3 = cimag(g->r1);   //a1
            r4 = cimag(g->r3);   //a2
            A = sqrt(sqr(r1-r2) + sqr(r3+r4));
            B = sqrt(sqr(r1-r2) + sqr(r3-r4));
            double g1 = sqrt((4.*sqr(r3)-sqr(A-B))/(sqr(A+B)-4.*sqr(r3)));
            mm = 4.*A*B/sqr(A+B);
            g->rp  = r1 - r3*g1;
            if (g->rp > r_bh(g->a)) warning("geodesic_priv_R_roots(): g->rp > r_bh(g->a)");
            g->Rpc = 2./(A+B) * jacobi_itn(-1./g1, mm);
            break;
        
        default:
            return FALSE;
    }


    return TRUE;
}




DEVICEFUNC
int geodesic_priv_T_roots(geodesic *g, double m, int *error)
// T integral roots
// M(m) = q + (a2 - l2 - q)m2 - a2*m4 = -a2*(m4 + m2*(q + l2 -a2)/a2 - q/a2) = a2(m2m+m2)(m2p-m2)
{
    double a  = g->a;
    double l  = g->l;
    double q  = g->q;
    double a2 = sqr(a);
    double l2 = sqr(l);

    // the straightforward solution of quadratic function fails for small spins due to numerical cancelation:
    //   g->m2m = 1./(2.*a2) * ( sqrt(sqr(qla)+4.*q*a2) + qla );
    //   g->m2p = 1./(2.*a2) * ( sqrt(sqr(qla)+4.*q*a2) - qla );
    // the roots can be calculated with the help of equality m2m*m2p=q/a^2
    // note: Dexter&Agol define m2m with opposite sign
    #ifndef CUDA
    long double qla = q + l2 - a2;
    long double X = sqrt(sqr(qla)+4.*q*a2) + qla;
    long double dbla = a2+a2;
    long double dblq = q+q;
    g->m2m = X/dbla;
    g->m2p = dblq/X;
    #else
    double qla = q + l2 - a2;
    double X = sqrt(sqr(qla)+4.*q*a2) + qla;
    double dbla = a2+a2;
    double dblq = q+q;
    g->m2m = X/dbla;
    g->m2p = dblq/X;
    #endif
    
    if ((g->m2p<=0.0) || (g->m2p >= 1.0)) {
        if (error) *error = GD_ERROR_MUPLUS_RANGE;
        return FALSE;
    }

    if (q > 0.0) {
        g->mm = g->m2p/(g->m2p+g->m2m);

        if ((g->mm<0.0) || (g->mm>=1.0)) {
            if (error) *error = GD_ERROR_MM_RANGE;
            return FALSE;
        }
        
        if (fabs(m) > sqrt(g->m2p)) {
            if (error) *error = GD_ERROR_MU0_RANGE;
            return FALSE;
        }

        g->mK = 1./sqrt(a2*(g->m2p+g->m2m));
    } else
    if (q < 0.0) {
        g->mm = (g->m2p+g->m2m)/g->m2p;

        if ((g->mm<0.0) || (g->mm>=1.0)) {
            if (error) *error = GD_ERROR_MM_RANGE;
            return FALSE;
        }

        if ((fabs(m) > sqrt(g->m2p)) || (fabs(m) < sqrt(-g->m2m))) {
            if (error) *error = GD_ERROR_MU0_RANGE;
            return FALSE;
        }

        g->mK = 1./sqrt(a2*g->m2p);
    } else {
        // q=0 means 
        //1. motion in the equatorial plane
        //2. motion off equatorial plane terminating in the singilarity
        if (error) *error = GD_ERROR_Q_RANGE;
        return FALSE;
    }

    return TRUE;
}


#undef theta_int
#undef theta_inv
