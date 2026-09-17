//************************************************************************
//    SIM5 library
//    sim5kerr-geod.c - null-geodesic motion routines
//------------------------------------------------------------------------
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (c) Michal Bursa, Astronomical Institute of the CAS
//************************************************************************


//! \file sim5kerr-geod.c
//! Null geodesics in Kerr spacetime.
//! 
//! Routines for computing null geodesics (trajectories of light rays) in Kerr spacetime based 
//! on numerical evaluation of elliptic integrals.


/*
#ifdef CUDA
__host__ __device__ void error(char *s) {
    //abort();
}
#endif
*/

//! \cond SKIP
// define some helper macros
#define theta_int(x) (g->mK*jacobi_icn((x)/sqrt(g->m2p),g->mm))
#define theta_inv(x) (sqrt(g->m2p)*jacobi_cn((x)/g->mK,g->mm))

// unit-private function declarations
DEVICEFUNC double geodesic_priv_RR(geodesic *g, double r);
DEVICEFUNC double geodesic_priv_TT(geodesic *g, double m);
DEVICEFUNC int geodesic_priv_R_roots(geodesic *g, double r0, int *error);
DEVICEFUNC int geodesic_priv_T_roots(geodesic *g, double m0, int *error);
//! \endcond



DEVICEFUNC
int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *error)
//! Initialization of a geodesic based on its impact parameters at infinity.
//! Makes a setup for a geodesic that is specified by its impact parameters at infinity. 
//! The imapct parameter is a perpendicular distance of the ray from the line-of-sight of the 
//! observer towards the black hole (the ray is parallel to the line-of-sight at infinity).
//!
//! @param i      inclination angle of observer (angle between BH rotation axis and direction to observer) [radians]
//! @param a      BH spin [0..1]
//! @param alpha  impact parameter in horizontal direction [GM/c^2]
//! @param beta   impact parameter in vertical direction [GM/c^2]
//! @param g      structure with information about geodesic (output)
//! @param error  error code (output)
//!
//! @result Returns TRUE on success or FALSE on error. In the latter case, an non-zero error code is 
//!         returned in `error` parameter. Information about the geodesic is stored in structure g.
{

    // validate the observer and spin before constructing trajectory constants.
    if ((a < 0.0) || (a > 1.-1e-6)) {
        if (error) *error = GD_ERROR_SPIN_RANGE;
        return FALSE;
    }

    if ((i <= 0.0) || (i>=PI_half)) {
        if (error) *error = GD_ERROR_INCL_RANGE;
        return FALSE;
    }

    // retain the existing regularization of a precisely horizontal image-plane ray.
    if (beta == 0.0) beta = +1e-6;

    g->a = fmax(1e-4, a);
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

    // vortical motion stays in one hemisphere and needs a dn phase, not cn.
    if (g->q<0.0) {
        double phase=geodesic_priv_dn_phase(g,g->cos_i);
        if (g->beta*g->cos_i>0.0) phase=-phase;
        g->Tip=phase*g->mK;
        g->Tpp=2.0*g->mK*elliptic_k(g->mm);
        if (error) *error=GD_OK;
        return TRUE;
    }

    // value of T-integral between turning points \int[-\mu_plus..\mu_plus]
	g->Tpp = 2.*theta_int(0.0);

    // value of T-integral between observer's position and turning point \int[cos_i..\mu_plus]
	g->Tip = theta_int(g->cos_i);

    if (error) *error = GD_OK;
    return TRUE;
}




DEVICEFUNC
int geodesic_init_src(double a, double r, double m, double k[4], int ppc, geodesic *g, int *error)
//! Initialization of a geodesic based on a position and direction.
//! Makes a setup for a geodesic that is specified by a point and a direction (4-momentum vector) 
//! somewhere along the trajectory. 
//!
//! @param a      BH spin [0..1]
//! @param r      radial coordinate of the point [GM/c^2]
//! @param m      poloidal coordinate of the point (\f$m=cos(\theta)\f$)
//! @param k      4-momentum null vector (\f$k \cdot k=0\f$) pointing in the direction of the ray
//! @param ppc    position with respect to pericenter (0=before pericenter, 1=after pericenter)
//! @param g      structure with information about geodesic (output)
//! @param error  error code (output)
//!
//! @result Returns TRUE on success or FALSE on error. In the latter case, an non-zero error code is 
//!         returned in `error` parameter. Information about the geodesic is stored in structure g.
{
    // reject sources for which the local coordinate description is not usable.
    if (!isfinite(a) || a<0.0 || a>=1.0 || !isfinite(r) ||
        r<=r_bh(a) || !isfinite(m) || fabs(m)>=1.0) {
        if (error) *error=GD_ERROR_ARGUMENT;
        return FALSE;
    }
    // calculate motion constants
    double l,q;
    photon_motion_constants(a, r, m, k, &l, &q);

    g->a = a;
    g->Q = 0.0;
    g->l = l;
    g->q = q;
    if (!isfinite(l) || !isfinite(q)) { if (error) *error=GD_ERROR_ARGUMENT; return FALSE; }

    // retain the emission momentum and handle the exact radial Schwarzschild limit.
    for (int j=0; j<4; ++j) g->k[j]=k[j];
    if (a==0.0 && l==0.0 && q==0.0) {
        g->type=GEOD_TYPE_RR_DBL; g->nrr=4;
        g->r1=g->r2=g->r3=g->r4=makeComplex(0.0,0.0);
        g->rp=0.0; g->Rpc=INFINITY;
        g->m2p=m*m; g->m2m=g->mm=0.0; g->mK=g->Tpp=INFINITY;
        g->p=1.0/r; g->cos_i=m; g->incl=acos(m);
        g->alpha=g->beta=g->Tip=0.0;
        if (error) *error=GD_OK;
        return TRUE;
    }
    // initialize observer quantities only after the radial and polar roots are valid.
    g->cos_i = g->alpha = g->beta = NAN;

    if (!geodesic_priv_R_roots(g, r, error)) return FALSE;
    if (!geodesic_priv_T_roots(g, m, error)) return FALSE;
    
    // store the source parameter using the appropriate outer or inner radial branch.
    g->p=geodesic_P_int(g,r,g->type==GEOD_TYPE_RR_BH ? k[1]<0.0 : ppc);
    if (!isfinite(g->p)) { if (error) *error=GD_ERROR_NUMERICAL; return FALSE; }
    // anchor the vortical phase at infinity while preserving the source direction.
    if (q<0.0) {
        double direction=(g->type==GEOD_TYPE_RC || g->type==GEOD_TYPE_CC) && k[1]>=0.0 ? -1.0 : 1.0;
        double dm_dP=(k[2]<0.0 ? 1.0 : -1.0)*direction;
        double phase=geodesic_priv_dn_phase(g,m);
        if (dm_dP*m>0.0) phase=-phase;
        // use an unwrapped phase so time integrals count complete polar oscillations.
        phase-=g->p/g->mK;
        g->Tip=phase*g->mK;
        g->Tpp=2.0*g->mK*elliptic_k(g->mm);
        // recover the observer coordinates without crossing to the opposite hemisphere.
        g->cos_i=copysign(sqrt(g->m2p)*geodesic_priv_dn(g,phase),m);
        g->incl=acos(g->cos_i);
        g->alpha=-l/sqrt(1.0-sqr(g->cos_i));
        double folded=phase-floor(phase/(2.0*elliptic_k(g->mm)))*2.0*elliptic_k(g->mm);
        double sign_dm=copysign(1.0,m)*(folded<elliptic_k(g->mm) ? -1.0 : 1.0);
        g->beta=sign_dm*sqrt(fmax(0.0,q-sqr(g->cos_i)*(sqr(g->alpha)-a*a)));
        if (error) *error=GD_OK;
        return TRUE;
    }

    // match the cn phase and its direction to the actual source, including outgoing RC rays.
    double direction=(g->type==GEOD_TYPE_RC || g->type==GEOD_TYPE_CC) && k[1]>=0.0 ? -1.0 : 1.0;
    double dm_dP=(k[2]<0.0 ? 1.0 : -1.0)*direction;
    double phase=jacobi_icn(fmax(-1.0,fmin(1.0,m/sqrt(g->m2p))),g->mm);
    if (dm_dP>0.0) phase=-phase;
    phase-=g->p/g->mK;
    // reduce the phase while retaining the sign that determines the observer's beta.
    double K=elliptic_k(g->mm);
    phase-=floor((phase+2.0*K)/(4.0*K))*4.0*K;
    g->Tpp=2.0*K*g->mK;
    g->Tip=fabs(phase)*g->mK;
    g->cos_i=sqrt(g->m2p)*jacobi_cn(fabs(phase),g->mm);
    g->incl=acos(g->cos_i);
    g->alpha=-l/sqrt(1.0-sqr(g->cos_i));
    g->beta=(phase<0.0 ? 1.0 : -1.0)*sqrt(fmax(0.0,q-sqr(g->cos_i)*(sqr(g->alpha)-a*a)));

    if (error) *error = GD_OK;
    return TRUE;
}




DEVICEFUNC
double geodesic_P_int(geodesic *g, double r, int ppc)
//! Position integral along geodesics at radius r.
//! It gives the value of the integral (Bursa 2017, eq. 34, and 43)
//! \f[ P = \int 1/\sqrt{R} dr = \int 1/\sqrt{\Theta} d\theta \f]
//! The integral is integrated from infinity to the given point on the trajecotry, 
//! where the geodesic reaches radius \f$r\f$ either before or behind the trajecory pericenter.
//!
//! The value of the integral increases monotonicly from infinity along the geodesic
//! and thus it provides a convenient way of parametrizing the position along the geodesic that is 
//! used in many other routined of the module. Note however, that the value of this integral is not 
//! the affine parameter, which would be another choice for parametrization.
//!
//! @param g    geodesic data
//! @param r    radial coordinate [GM/c^2]
//! @param ppc  position with respect to pericenter (0=before pericenter, 1=after pericenter)
//!
//! @result     Value of the position integral between infinity and given point.
{
    // the exact radial limit has a simple positional integral with no turning point.
    if (g->a==0.0 && g->l==0.0 && g->q==0.0) return 1.0/r;
    // compactify the CC integral at infinity instead of inventing a pericenter.
    if (g->type==GEOD_TYPE_CC) {
        if (!(r>r_bh(g->a))) return NAN;
        return geodesic_priv_cc_integral(g,0.0,1.0/r,0);
    }

    double r1,r2,r3,r4,u,v, mm, R, A, B;

    // the inner four-root branch lies below its upper turning radius by construction.
    #ifdef CUDA
    if (r < g->rp && g->type!=GEOD_TYPE_RR_BH) asm("exit;");
    #else
    if (r < g->rp && g->type!=GEOD_TYPE_RR_BH) error("(geodesic_P_int): r < periastron (%.3e < %.3e; nrr=%d)", r, g->rp, g->nrr);
    #endif
    if (r == g->rp) return g->Rpc;


    switch (g->type) {
        case GEOD_TYPE_RR:
            // distinguish the ingoing and outgoing portions of the outer allowed region.
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt(((r2-r4)*(r-r1))/((r1-r4)*(r-r2))), mm);
            return (ppc) ? g->Rpc + R : g->Rpc - R;
        
        case GEOD_TYPE_RR_DBL:
            #ifndef CUDA
            error("(geodesic_P_int): not implemented for GEOD_TYPE_RR_DBL");
            #endif
            return NAN;

        case GEOD_TYPE_RR_BH:
            // measure the inner branch relative to its upper radial turning point.
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            R  = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r1-r3)/(r2-r3)*(r2-r)/(r1-r)), mm);
            return (ppc) ? g->Rpc + R : g->Rpc - R;

        case GEOD_TYPE_RC:
            // the exterior RC branch is monotonic, so its parameter runs from infinity inward.
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            u  = creal(g->r3);
            v  = cimag(g->r3);
            A = sqrt(sqr(r1-u)+sqr(v));
            B = sqrt(sqr(r2-u)+sqr(v));
            mm = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
            R  = 1./sqrt(A*B) * jacobi_icn(((A-B)*r+r1*B-r2*A)/((A+B)*r-r1*B-r2*A), mm);
            //if (R > g->Rpc) fprintf(stderr, "RC: res=%e Rpc=%e R=%e r=%.4e r1=%.4e r2=%.4e A=%.4e B=%.4e mm=%.4e z=%.4e icn=%.4e\n", g->Rpc-R, g->Rpc, R, r, r1, r2, A, B, mm, ((A-B)*r+r1*B-r2*A)/((A+B)*r-r1*B-r2*A), jacobi_icn(((A-B)*r+r1*B-r2*A)/((A+B)*r-r1*B-r2*A), mm));
            // the RC case has no exterior turning point, so 0 < P < Rpc.
            return g->Rpc-R;


    }

    // this point should never be reached
    return NAN;
}



//! \cond SKIP
DEVICEFUNC
void geodesic_position(geodesic *g, double P, double x[4])
//! Point on the geodesic, where the position integral gains value P.
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
//! \endcond



DEVICEFUNC
double geodesic_position_rad(geodesic *g, double P)
//! Radius at which the position integral gains value P.
//! Given the value \f$P\f$ of the positional integral along the geodesic, the function computes
//! the radial coordinate for the position.
//!
//! @param g  geodesic data
//! @param P  value of the position integral
//!
//! @result Radial coordinate value [GM/c^2] or NAN in case of error.
{
    // invert the exact radial limit without introducing elliptic singularities.
    if (g->a==0.0 && g->l==0.0 && g->q==0.0) return P>0.0 ? 1.0/P : NAN;
    // invert the monotonic exterior CC integral on a finite inverse-radius bracket.
    if (g->type==GEOD_TYPE_CC) {
        double lo=0.0,hi=1.0/r_bh(g->a);
        double maxP=geodesic_priv_cc_integral(g,0.0,hi,0);
        if (!(P>0.0 && P<maxP)) return NAN;
        // retain the bracket until the radius is resolved to floating-point precision.
        for (int i=0; i<60; ++i) {
            double u=(lo+hi)/2.0;
            double value=geodesic_priv_cc_integral(g,0.0,u,0);
            if (!isfinite(value)) return NAN;
            if (value<P) lo=u; else hi=u;
        }
        return 2.0/(lo+hi);
    }

    double r1,r2,r3,r4,u,v;

    // reject parameters outside the established analytic branch domain.
    if ((P<=0.0)||(P>=2.*g->Rpc)) {
        #ifndef CUDA
        error("(geodesic_position_rad) P out of range (P=%e, 2Rpa=%e)", P, 2*g->Rpc);
        #endif
        return NAN;
    }
    if (P == g->Rpc) return g->rp;

    switch (g->type) {
        case GEOD_TYPE_RR:
            // the absolute phase covers both sides of the outer radial turning point.
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            double m4 = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            double x4 = 0.5*fabs    (P - g->Rpc)*sqrt((r2-r4)*(r1-r3));
            double sn2 = pow( jacobi_sn(x4,m4), 2.0);
            return ( r1*(r2-r4)-r2*(r1-r4)*sn2 ) / ( r2-r4-(r1-r4)*sn2 );

        case GEOD_TYPE_RR_DBL:
            #ifndef CUDA
            error("(geodesic_position_rad): not implemented for GEOD_TYPE_RR_DBL");
            #endif
            return NAN;

        case GEOD_TYPE_RR_BH:
            #ifndef CUDA
            error("(geodesic_position_rad): not implemented for GEOD_TYPE_RR_BH");
            #endif
            return NAN;

        case GEOD_TYPE_RC:
            // invert the monotonic mixed-root branch using its Jacobi cn solution.
            // the RC case has no exterior turning point, so P cannot exceed Rpc.
            if (P > g->Rpc) return NAN;
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            u  = creal(g->r3);
            v  = cimag(g->r3);
            double A = sqrt(sqr(r1-u)+sqr(v));
            double B = sqrt(sqr(r2-u)+sqr(v));
            double m2 = (sqr(A+B) - sqr(r1-r2)) / (4.*A*B);
            double cn = jacobi_cn(sqrt(A*B)*(g->Rpc-P), m2);
            //fprintf(stderr, "r-RC: mm=%.4e z=%.4e icn=%.4e  r=%.4e\n", m2, sqrt(A*B)*(g->Rpc-P), cn, (r2*A - r1*B - (r2*A+r1*B)*cn ) / ( (A-B) - (A+B)*cn ));
            return (r2*A - r1*B - (r2*A+r1*B)*cn ) / ( (A-B) - (A+B)*cn );
        

    }

    // this point should never be reached
    return NAN;
}




DEVICEFUNC
double geodesic_position_pol(geodesic *g, double P)
//! Poloidal coordinate value at which the position integral gains value P.
//! Given the value \f$P\f$ of the positional integral along the geodesic, the function computes
//! the poloidal coordinate for the position. The coordinate is returned as a cosine of the angle theta.
//!
//! @param g  geodesic data
//! @param P  value of the position integral
//!
//! @result Poloidal coordinate value [cos(theta)] or NAN in case of error.
{
    // vortical rays retain their hemisphere while the dn phase repeats.
    if (g->q<0.0)
        return copysign(sqrt(g->m2p)*geodesic_priv_dn(g,(g->Tip+P)/g->mK),g->cos_i);
    // a radial Schwarzschild photon keeps its initial polar angle.
    if (g->a==0.0 && g->l==0.0 && g->q==0.0) return g->cos_i;

    double sign_dm, T;

    switch (g->type) {
        // trajectories that go to infinity
        case GEOD_TYPE_RR:
        case GEOD_TYPE_RC:
        case GEOD_TYPE_CC:
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

        case GEOD_TYPE_RR_DBL:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for GEOD_TYPE_RR_DBL");
            #endif
            return NAN;

        case GEOD_TYPE_RR_BH:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for GEOD_TYPE_BH");
            #endif
            return NAN;
    }

    // default
    return NAN;
}




DEVICEFUNC
double geodesic_position_pol_sign_k_theta(geodesic *g, double P)
//! Sign of the \f$k^\theta\f$ component of the 4-momentum.
//! Gives the orientation of the 4-momentum vector in the poloidal direction by 
//! returning the sign of \f$k^\theta\f$ component of the momentum vector.
//!
//! @param g  geodesic data
//! @param P  value of the position integral
//!
//! @result Returns +1 or -1 or NAN in case of an error.
{
    double sign_dm, T;

    switch (g->type) {
        // trajectories that go to infinity
        case GEOD_TYPE_RR:
        case GEOD_TYPE_RC:
        case GEOD_TYPE_CC:
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

        case GEOD_TYPE_RR_DBL:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for GEOD_TYPE_RR_DBL");
            #endif
            return NAN;

        case GEOD_TYPE_RR_BH:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for GEOD_TYPE_BH");
            #endif
            return NAN;
    }

    // default
    return NAN;
}




DEVICEFUNC static double geodesic_priv_azm_phase(double phase, double complement, double mm)
//! Unwrapped third-kind elliptic primitive for the polar azimuth integral.
//! - integrates 1/(1-n*sn(x|mm)^2) from x=0 to an arbitrary signed phase
//! - counts full 2K periods explicitly so complete azimuthal windings are retained
//! - uses Carlson integrals with an explicit 1-n to retain precision near the axis
//! @param phase dimensionless Jacobi phase, including any complete polar periods
//! @param complement positive 1-n, where the elliptic characteristic n is in [0,1)
//! @param mm elliptic parameter in [0,1)
//! @result Signed dimensionless primitive; the caller supplies the angular scale.
{
    // retain complete periods before evaluating the principal elliptic integral.
    double K=elliptic_k(mm), n=1.0-complement;
    double complete=rf(0.0,1.0-mm,1.0)+n*rj(0.0,1.0-mm,1.0,complement)/3.0;
    double cycles=floor(phase/(2.0*K)), rem=phase-cycles*2.0*K;
    // keep cn squared and 1-n explicit to avoid subtracting nearly equal numbers.
    double sn,cn,dn;
    jacobi_sncndn(fmin(rem,2.0*K-rem),mm,&sn,&cn,&dn);
    double s2=sn*sn, c2=cn*cn;
    double partial=sn*(rf(c2,dn*dn,1.0)+n*s2*rj(c2,dn*dn,1.0,c2+complement*s2)/3.0);
    if (rem==K) partial=complete;
    // reflection about K handles either half of a polar period continuously.
    return 2.0*cycles*complete+(rem>K ? 2.0*complete-partial : partial);
}

DEVICEFUNC
double geodesic_position_azm(geodesic *g, double r, double m, double P)
//! Unwrapped azimuthal primitive measured from infinity (P=0).
//! - integrates dphi/dP=(2*a*r-l*a^2)/Delta+l/(1-m^2), Bursa (2017), Eq. (20)
//! - includes radial turning points for RR and polar turning points for both signs of q
//! - uses elliptic radial integrals for RR/RC and convergent radial quadrature for CC
//! - supports zero spin, including the exact radial Schwarzschild limit
//! @param g initialized geodesic with its polar phase anchored at P=0
//! @param r exterior Boyer-Lindquist radius; zero requests inference of both r and m
//! @param m cos(theta) at P; must correspond to the same trajectory position
//! @param P nonnegative positional integral, with P=0 denoting the reference infinity
//! @result Signed angle in radians, without reduction modulo 2*pi; NaN for invalid
//! or unsupported positions. For future motion with decreasing P (outgoing RC/CC),
//! negate the difference of endpoint primitives to obtain the physical azimuth change.
//! Azimuth is a singular coordinate on the axis; an exactly axial crossing has no
//! uniquely defined continuous azimuth there.
{
    // allow the reference infinity and infer coordinates only when explicitly omitted.
    if (!isfinite(P) || P<0.0) return NAN;
    if (P==0.0) return 0.0;
    if (r==0.0) {
        r=geodesic_position_rad(g,P);
        m=geodesic_position_pol(g,P);
    }
    if (!isfinite(r) || r<=r_bh(g->a) || !isfinite(m) || fabs(m)>=1.0) return NAN;
    // zero spin has no radial frame dragging, avoiding singular horizon primitives.
    double phi=0.0, a2=sqr(g->a);
    if (g->a!=0.0) {
        double rp=1.0+sqrt(1.0-a2), rm=1.0-sqrt(1.0-a2);
        double r1=creal(g->r1), r2=creal(g->r2), A, B;
        switch (g->type) {
            case GEOD_TYPE_RR: {
                // reverse the radial leg after pericenter while keeping P increasing.
                double sign=P>g->Rpc ? 1.0 : -1.0;
                double r3=creal(g->r3), r4=creal(g->r4);
                A=integral_R_rp_re_inf(r1,r2,r3,r4,rp)+sign*integral_R_rp_re(r1,r2,r3,r4,rp,r);
                B=integral_R_rp_re_inf(r1,r2,r3,r4,rm)+sign*integral_R_rp_re(r1,r2,r3,r4,rm,r);
                phi=(A*(g->a*rp-g->l*a2/2.0)-B*(g->a*rm-g->l*a2/2.0))/sqrt(1.0-a2);
                break;
            }
            case GEOD_TYPE_RC:
                // the exterior RC branch is monotonic between infinity and the horizon.
                A=integral_R_rp_cc2_inf(r1,r2,g->r3,rp,r);
                B=integral_R_rp_cc2_inf(r1,r2,g->r3,rm,r);
                phi=(A*(g->a*rp-g->l*a2/2.0)-B*(g->a*rm-g->l*a2/2.0))/sqrt(1.0-a2);
                break;
            case GEOD_TYPE_CC:
                // inverse radius makes the CC integral from infinity finite and smooth.
                phi=geodesic_priv_cc_integral(g,0.0,1.0/r,2);
                break;
            default:
                return NAN;
        }
    }
    // zero angular momentum has no polar azimuth rate, even at zero spin.
    if (g->l==0.0) return phi;
    // recover 1-m2p from the polar potential at the axis, without cancellation.
    double complement=g->a==0.0 ? sqr(g->l)/(g->q+sqr(g->l))
                               : sqr(g->l)/(a2*(1.0+g->m2m));
    if (!(complement>0.0)) return NAN;
    // shift cn/dn by K so the third-kind primitive has a positive characteristic.
    double phase0=(g->q<0.0 ? g->Tip : (g->beta>=0.0 ? -g->Tip : g->Tip))/g->mK;
    phase0+=elliptic_k(g->mm);
    double phase=phase0+P/g->mK;
    double scale, linear, pi_complement;
    if (g->q<0.0) {
        double d=complement+g->m2p*g->mm;
        linear=1.0;
        scale=g->m2p*(1.0-g->mm)/d;
        pi_complement=complement*(1.0-g->mm)/d;
    } else {
        double n=g->mm+g->m2p*(1.0-g->mm);
        linear=g->mm/n;
        scale=g->m2p*(1.0-g->mm)/n;
        pi_complement=complement*(1.0-g->mm);
    }
    // add the unwrapped polar contribution without losing complete turns.
    phi+=g->l*(linear*P+g->mK*scale*(geodesic_priv_azm_phase(phase,pi_complement,g->mm)
                                  -geodesic_priv_azm_phase(phase0,pi_complement,g->mm)));
    return phi;
}



/* V-phase for the polar time integral: int from m up to +sqrt(m2p) of
 * m^2 dm/sqrt((m2p-m^2)(m^2+m2m)).  Rises as m falls, exactly as the polar
 * phase theta_int does, so the two share a bounce structure and one fold
 * serves both.  m^2 is even, hence the reflection for m < 0. */
DEVICEFUNC
double geodesic_priv_vphase(geodesic *g, double m)
{
    double J0 = integral_T_m2(g->m2m, g->m2p, 0.0);
    double Jx = integral_T_m2(g->m2m, g->m2p, fabs(m));
    return (m >= 0.0) ? Jx : 2.0*J0 - Jx;
}


DEVICEFUNC
double geodesic_timedelay(geodesic *g, double P1, double r1, double m1, double P2, double r2, double m2)
//! Time delay (travel time) between positions P1 and P2.
//! Gives time it takes the light to travel between two points along a geodesic.
//! Returned value is always positive independent of the relative position of 
//! P1 and P2 along the geodesic.
//!
//! @param g   geodesic data
//! @param P1  value of the position integral at point A
//! @param r1  value of the radial coordinate at point A; if zero, it is computed from P1 internally
//! @param m1  cos(theta) at point A; inferred together with r1 only when r1=0
//! @param P2  value of the position integral at point B
//! @param r2  value of the radial coordinate at point B; if zero, it is computed from P2 internally
//! @param m2  cos(theta) at point B; inferred together with r2 only when r2=0
//!
//! @result Nonnegative coordinate time in GM/c^3, or NaN for invalid/unsupported endpoints.
//! Endpoints must lie outside the horizon. RR/RC use elliptic radial integrals;
//! CC uses convergent radial quadrature. Vortical (q<0) polar motion is included.
{
    double time = 0.0;

    // exchange endpoints if needed, so it can be assumed P1 < P2
    if (P1 > P2) {
        double tmp;
        tmp=P2; P2=P1; P1=tmp;
        tmp=r2; r2=r1; r1=tmp;
        tmp=m2; m2=m1; m1=tmp;
    }

    // infer both coordinates only when the caller omits the radius.
    if (r1 == 0) {
        r1 = geodesic_position_rad(g, P1);
        m1 = geodesic_position_pol(g, P1);
    }

    if (r2 == 0) {
        r2 = geodesic_position_rad(g, P2);
        m2 = geodesic_position_pol(g, P2);
    }


    // coordinate time is defined here only between finite exterior endpoints.
    if (!isfinite(P1) || !isfinite(P2) || !isfinite(r1) || !isfinite(r2) ||
        !isfinite(m1) || !isfinite(m2) || r1<=r_bh(g->a) || r2<=r_bh(g->a)) return NAN;
    if (P1==P2) return 0.0;
    // use the exact radial Schwarzschild delay rather than singular elliptic limits.
    if (g->a==0.0 && g->l==0.0 && g->q==0.0)
        return fabs((r2-r1)+2.0*log((r2-2.0)/(r1-2.0)));

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

    switch (g->type) {
        // trajectories that go to infinity
        case GEOD_TYPE_RR:
            // add radial legs across a turning point, otherwise subtract their primitives.
            s = (((P1 > g->Rpc)&&(P2 < g->Rpc)) || ((P1 < g->Rpc)&&(P2 > g->Rpc))) ? +1 : -1;
            R0 = integral_R_r0_re(ra, rb, rc, rd, r1)     + s*integral_R_r0_re(ra, rb, rc, rd, r2);
            R1 = integral_R_r1_re(ra, rb, rc, rd, r1)     + s*integral_R_r1_re(ra, rb, rc, rd, r2);
            R2 = integral_R_r2_re(ra, rb, rc, rd, r1)     + s*integral_R_r2_re(ra, rb, rc, rd, r2);
            RA = integral_R_rp_re(ra, rb, rc, rd, rp, r1) + s*integral_R_rp_re(ra, rb, rc, rd, rp, r2);
            // the inner-horizon term vanishes at zero spin; do not evaluate its singular primitive.
            RB = a2==0.0 ? 0.0 : integral_R_rp_re(ra, rb, rc, rd, rm, r1) + s*integral_R_rp_re(ra, rb, rc, rd, rm, r2);
            A = (-g->a*g->l+4.)*rp - 2.*a2;
            B = (+g->a*g->l-4.)*rm + 2.*a2;
            time += 4.*fabs(R0) + 2.*fabs(R1) + fabs(R2) + (A*fabs(RA) + B*fabs(RB))/sqrt(1.-a2);
            //fprintf(stderr,"RR1=%e  RR2=%e\n", geodesic_priv_RR(g,r1), geodesic_priv_RR(g,r2));
            //fprintf(stderr,"ra=%.3e rb=%.3e rc=%.3e rd=%.3e\n", ra, rb, rc, rd);
            //fprintf(stderr,"R0=%.3e R1=%.3e R2=%.3e RA=%.3e RB=%.3e A=%.3e B=%.3e dt=%.3e \n", R0,R1,R2,RA,RB,A,B,time);
            break;

        case GEOD_TYPE_RR_DBL:
            #ifndef CUDA
            error("(geodesic_timedelay): not implemented for GEOD_TYPE_RR_DBL");
            #endif
            return NAN;

        case GEOD_TYPE_RR_BH:
            #ifndef CUDA
            error("(geodesic_timedelay): not implemented for GEOD_TYPE_RR_BH");
            #endif
            return NAN;

        case GEOD_TYPE_RC:
            // this exterior radial branch is monotonic, so use endpoint differences.
            R0 = integral_R_r0_cc(ra, rb, g->r3, r1) -integral_R_r0_cc(ra, rb, g->r3, r2);
            R1 = (r1<r2) ? integral_R_r1_cc(ra, rb, g->r3, r1, r2) : integral_R_r1_cc(ra, rb, g->r3, r2, r1);
            R2 = (r1<r2) ? integral_R_r2_cc(ra, rb, g->r3, r1, r2) : integral_R_r2_cc(ra, rb, g->r3, r2, r1);
            RA = (r1<r2) ? integral_R_rp_cc2(ra, rb, g->r3, rp, r1, r2) : integral_R_rp_cc2(ra, rb, g->r3, rp, r2, r1);
            // the inner-horizon term vanishes at zero spin; do not evaluate its singular primitive.
            RB = a2==0.0 ? 0.0 : (r1<r2) ? integral_R_rp_cc2(ra, rb, g->r3, rm, r1, r2) : integral_R_rp_cc2(ra, rb, g->r3, rm, r2, r1);
            A = (-g->a*g->l+4.)*rp - 2.*a2;
            B = (+g->a*g->l-4.)*rm + 2.*a2;
            time += 4.*fabs(R0) + 2.*fabs(R1) + fabs(R2) + (A*fabs(RA) + B*fabs(RB))/sqrt(1.-a2);
            //fprintf(stderr,"RR1=%e  RR2=%e\n", geodesic_priv_RR(g,r1), geodesic_priv_RR(g,r2));
            //fprintf(stderr,"ra=%.3e rb=%.3e rc=%.3e rd=%.3e\n", ra, rb, rc, rd);
            //fprintf(stderr,"R0=%.3e R1=%.3e R2=%.3e RA=%.3e RB=%.3e A=%.3e B=%.3e dt=%.3e \n", R0,R1,R2,RA,RB,A,B,time);
            break;

        case GEOD_TYPE_CC:
            // no radial turning point. Integrate Eq. (19) in log(r-r_h)
            // to resolve both near-horizon and large-radius endpoints.
            // the polar contribution is added below, including q<0.
            time=fabs(geodesic_priv_cc_integral(g,log(r1-rp),log(r2-rp),1));
            if (!isfinite(time)) return NAN;
            break;

    }

    // add the polar part of coordinate time; the radial primitives do not contain it.
    // fold the q>0 trajectory through each polar bounce and infer its direction from
    // the supplied endpoints, retaining the existing RR/RC phase convention.
    if ((g->a != 0.0) && (g->q > 0.0) && (g->Tpp > 0.0)) {
        double Tpp  = g->Tpp;
        double Vpp  = 2.0*integral_T_m2(g->m2m, g->m2p, 0.0);
        double P    = fabs(P2-P1);
        double tau1 = theta_int(m1);
        double psi1 = geodesic_priv_vphase(g, m1);
        double best = 0.0, best_miss = -1.0;
        int    s;

        for (s=0; s<2; s++) {
            double dir = s ? -1.0 : +1.0;
            double xi  = tau1 + dir*P;
            double n   = floor(xi/Tpp);
            double r   = xi - n*Tpp;
            long   odd = (((long)n) % 2 + 2) % 2;
            double tf  = odd ? (Tpp - r) : r;
            double me  = theta_inv(tf);
            double V   = geodesic_priv_vphase(g, me);
            double psi2 = n*Vpp + (odd ? (Vpp - V) : V);
            double miss = fabs(me - m2);
            if ((best_miss < 0.0) || (miss < best_miss)) {
                best_miss = miss;
                best = g->a * fabs(psi2 - psi1);
            }
        }
        if (isfinite(best)) time += best;
    }

    // preserve the actual dn phase for vortical rays, including complete oscillations.
    if (g->q<0.0) time+=geodesic_priv_vortical_time(g,P1,P2);
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

    switch (g->type) {
        // trajectories that go to infinity
        case GEOD_TYPE_RR:
        case GEOD_TYPE_RC:
        case GEOD_TYPE_CC:
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

        case GEOD_TYPE_RR_DBL:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for GEOD_TYPE_RR_DBL");
            #endif
            return NAN;

        case GEOD_TYPE_RR_BH:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for GEOD_TYPE_BH");
            #endif
            return NAN;
    }

    // default
    return NAN;
}




DEVICEFUNC
void geodesic_momentum(geodesic *g, double P, double r, double m, double k[])
//! Photon 4-momentum.
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
//! @param g      - geodesic data
//! @param P      - value of the position integral
//! @param r      - radial coordinate (value or zero)
//! @param m      - cosine of poloidal coordinate (value or zero)
//! @param k      - 4-momentum vector (output)
//!
//! @result Photon 4-momentum vector is returned in k[].
{
    double dm;

    // calc [r,m] coordinate if it has not been provided
    if ((r==0.0) && (m==0.0)) {
        r = geodesic_position_rad(g, P);
        m = geodesic_position_pol(g, P);
    }

    switch (g->type) {
        // trajectories that go to infinity
        case GEOD_TYPE_RR:
        case GEOD_TYPE_RC:
        case GEOD_TYPE_CC:
            dm = geodesic_dm_sign(g, P);
            photon_momentum(g->a, 0.0, r, m, g->l, g->q, (P<g->Rpc?-1:+1), dm, k);
            return;
            
        case GEOD_TYPE_RR_DBL:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for GEOD_TYPE_RR_DBL");
            #endif
            k[0]=k[1]=k[2]=k[3]=NAN;
            return;

        case GEOD_TYPE_RR_BH:
            #ifndef CUDA
            error("(geodesic_dm_sign): not implemented for GEOD_TYPE_BH");
            #endif
            k[0]=k[1]=k[2]=k[3]=NAN;
            return;
    }

    return;
}




DEVICEFUNC
double geodesic_find_midplane_crossing(geodesic *g, int order)
//! Finds a crossing of the geodesic with the equatorial plane.
//! Calculates, where (if ever) the geodesic crosses the equatorial plane, and returns
//! the value of the positional integral for that place. This is the fastest way to 
//! integrate over the equatorial plane. The positional integral can be converted to radius, which 
//! allows straightforward integration over the solid angle, for example
//! \f[ F_\nu(E) = 1/D^2 \int I_\nu(E/g, r) d\alpha\,d\beta \f]
//! where \f$r = r(\alpha, \beta)\f$, \f$D\f$ is distance and \f$g\f$ is g-factor.
//!
//! @param g      geodesic data
//! @param order  order of crossing; order=0 zero is the first crossing, higher orders may be reached 
//!                by some geodesics that loop around the photon orbit
//!
//! @result Value of the positional integral at the equatorial plane.
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
//! Makes a step along the geodesic.
//! Moves current position on the ray along the geodescis. On input, the function receives the 
//! current position integral value and radial and poloidal coordinate. It computes new values for 
//! P, r and m and shifts them to a new position along the geodesics by a given step.
//! 
//! This function is meant to be called in a cycle to follow the geodesics in a piecewise steps.
//!
//! @param g        geodesic data
//! @param step     size of step to advance
//! @param P        value of the positional integral (input and output)
//! @param r        value of the radial coordinate (input and output)
//! @param m        value of the poloidal coordinate (input and output; \f$m=cos(theta)\f$)
//! @param status   status code; status=0 if ok, it get a non-zero value on an error
{
    const double MAXSTEP_FACTOR = 5e-2;

    do {
        double truestep = step/fabs(step) * fmin(fabs(step), MAXSTEP_FACTOR*sqrt(*r));
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

//! \cond SKIP

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
//! Initialize the radial root classification and positional-integral coefficients.
//! @param g geodesic with a,l,q already initialized; roots and branch data are outputs
//! @param r0 source radius used to distinguish the outer and inner four-root regions
//! @param error optional output error code on failure
//! @result TRUE for a supported radial branch, FALSE for invalid or degenerate roots.
//! Real roots are bracketed before the remaining complex roots are reconstructed.
{
    double a  = g->a;
    double l  = g->l;
    double q  = g->q;
    double a2 = sqr(a);
    double l2 = sqr(l);
    // bracket real roots before reconstructing any complex pair.
    geod_equations eq={a,l,q,a2-l2-q,q+sqr(l-a),a2*q};
    double real[4];
    g->nrr=geodesic_priv_radial_roots(&eq,real);
    if (g->nrr==4) {
        g->r1=makeComplex(real[3],0.0); g->r2=makeComplex(real[2],0.0);
        g->r3=makeComplex(real[1],0.0); g->r4=makeComplex(real[0],0.0);
    } else if (g->nrr==2) {
        // recover the conjugate pair from coefficient identities after fixing the real roots.
        g->r1=makeComplex(real[1],0.0); g->r2=makeComplex(real[0],0.0);
        double u=-(real[0]+real[1])/2.0;
        double v2=eq.A-real[0]*real[1]+3.0*u*u;
        if (!(v2>0.0)) { if (error) *error=GD_ERROR_NUMERICAL; return FALSE; }
        g->r3=makeComplex(u,sqrt(v2)); g->r4=conj(g->r3);
    } else if (g->nrr==0) {
        // with no real roots there is no large/small real-root cancellation.
        int nr;
        quartic_eq_c(0.0,eq.A,2.0*eq.B,-eq.C,&nr,&g->r1,&g->r2,&g->r3,&g->r4);
        if (nr!=0) { if (error) *error=GD_ERROR_NUMERICAL; return FALSE; }
    } else {
        if (error) *error=GD_ERROR_TYPE_RR_DOUBLE;
        return FALSE;
    }

    // classify the allowed radial region containing the source.
    switch (g->nrr) {
        case 4:
            g->type = GEOD_TYPE_RR;
            // r0 can only be between r3 and r2; or it can be above r1
            // anything else is an error
            if ((r0<creal(g->r3)) || ((r0>creal(g->r2)) && (r0<creal(g->r1)))) {
                if (error) *error = GD_ERROR_UNKNOWN_SOLUTION;
                return FALSE;
            }
            // if r1 and r2 are close to each other, it is a double root solution
            if (fabs(creal(g->r1)-creal(g->r2)) < 1e-8) {
                g->type = GEOD_TYPE_RR_DBL;
                if (error) *error = GD_ERROR_TYPE_RR_DOUBLE;
                return FALSE;
            }
            // if r0 is between r3 and r2, it is an inner solution
            if ((r0>=creal(g->r3)) && (r0<=creal(g->r2))) {
                g->type = GEOD_TYPE_RR_BH;
            }
            break;
        case 2:
            g->type = GEOD_TYPE_RC;
            break;
        case 0:
            g->type = GEOD_TYPE_CC;
            break;
        default:
            if (error) *error = GD_ERROR_UNKNOWN_SOLUTION;
            return FALSE;
    }

    // prepare the branch constants needed by the existing positional integrals.
    double r1,r2,r3,r4,u,v,mm,A,B;
    switch (g->type) {
        case GEOD_TYPE_RR:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            g->rp  = r1;
            g->Rpc = 2./sqrt((r1-r3)*(r2-r4)) * jacobi_isn(sqrt((r2-r4)/(r1-r4)), mm);
            break;

        case GEOD_TYPE_RR_BH:
            r1 = creal(g->r1);
            r2 = creal(g->r2);
            r3 = creal(g->r3);
            r4 = creal(g->r4);
            mm = ((r2-r3)*(r1-r4))/((r2-r4)*(r1-r3));
            g->rp  = r2;
            g->Rpc = 2./sqrt((r1-r3)*(r2-r4)) * elliptic_k(mm); //ellipticK(mm) = jacobi_isn(1, mm)
            break;

        case GEOD_TYPE_RC:
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

        case GEOD_TYPE_CC:
            // no real root means no pericenter or finite turning-point integral.
            g->rp=-INFINITY;
            g->Rpc=INFINITY;
            break;

        default:
            return FALSE;
    }


    return TRUE;
}




DEVICEFUNC
int geodesic_priv_T_roots(geodesic *g, double m, int *error)
//! Initialize polar turning points and the Jacobi phase scaling.
//! @param g geodesic with a,l,q initialized; root and elliptic parameters are outputs
//! @param m source or observer cos(theta), required to lie in the allowed polar range
//! @param error optional output error code on failure
//! @result TRUE when the polar phase is supported, FALSE for an invalid or degenerate
//! polar potential. The Schwarzschild q>0 limit is handled without dividing by a^2.
{
    double a  = g->a;
    double l  = g->l;
    double q  = g->q;
    double a2 = sqr(a);
    double l2 = sqr(l);

    // take the Schwarzschild limit explicitly to avoid division by a squared.
    if (a==0.0) {
        if (!(q>0.0)) { if (error) *error=GD_ERROR_Q_RANGE; return FALSE; }
        g->m2p=q/(q+l2); g->m2m=INFINITY; g->mm=0.0;
        g->mK=1.0/sqrt(q+l2);
        if (m*m>g->m2p+1e-12) { if (error) *error=GD_ERROR_MU0_RANGE; return FALSE; }
        return TRUE;
    }

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
    
    // allow the axis as a polar turning point, including roundoff near m squared equals one.
    if ((g->m2p<=0.0) || (g->m2p > 1.0+1e-12)) {
        if (error) *error = GD_ERROR_MUPLUS_RANGE;
        return FALSE;
    }

    // ordinary polar motion can pass through the equatorial plane.
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
    // vortical motion stays between two turning points in the same hemisphere.
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

//! \endcond


#undef theta_int
#undef theta_inv
