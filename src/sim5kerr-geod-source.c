//! \file sim5kerr-geod-source.c
//! Source intersections and numerical support for Kerr geodesics.
//! Private helpers are shared with sim5kerr-geod.c in the amalgamated library.

//! \cond SKIP
enum { GEOD_U, GEOD_V, GEOD_M, GEOD_W, GEOD_N };

typedef struct {
    double a, l, q, A, B, C;
} geod_equations;

// bursa (2017), Eqs. (21)-(22), with u=1/r. Energy-rescaled Mino time
// (not the paper's affine tau) removes coordinate singularities from
// radial/polar motion. Differentiating the squared first integrals gives
// smooth second-order equations, including at turning points:
// u'^2 = 1 + A*u^2 + 2*B*u^3 - C*u^4
// m'^2 = q + A*m^2 - a^2*m^4
// A=a^2-l^2-q, B=q+(l-a)^2, C=a^2*q.
DEVICEFUNC static void geodesic_priv_rhs(const geod_equations *p, const double y[GEOD_N], double d[GEOD_N])
//! Derivatives of the separated radial and polar motion in energy-rescaled Mino time.
//! - uses u=1/r and m=cos(theta), with their first derivatives v and w
//! - evolves the differentiated potentials smoothly through turning points
//! - does not evolve coordinate time; geodesic_timedelay() owns that calculation
//! @param p immutable Kerr constants and radial/polar polynomial coefficients
//! @param y state ordered as (u,v,m,w); lengths are in GM/c^2
//! @param d output derivatives in the same order
//! @result Fills d without changing the equations or input state.
{
    // evolve radial motion without choosing square-root signs at turning points.
    double u = y[GEOD_U], m = y[GEOD_M], a2 = p->a*p->a;
    d[GEOD_U] = y[GEOD_V];
    d[GEOD_V] = u*(p->A + u*(3.0*p->B - 2.0*p->C*u));
    // the polar acceleration preserves the separated angular potential.
    d[GEOD_M] = y[GEOD_W];
    d[GEOD_W] = m*(p->A - 2.0*a2*m*m);
    // geometry only: coordinate time belongs to geodesic_timedelay().

}

DEVICEFUNC static void geodesic_priv_rk4(const geod_equations *p, const double y[GEOD_N], double h, double out[GEOD_N])
//! One fourth-order Runge-Kutta step for separated geodesic motion.
//! - evaluates four stages without changing the supplied initial state
//! - step acceptance and error control are the responsibility of the caller
//! @param p immutable constants defining the geodesic
//! @param y initial (u,v,m,w) state
//! @param h step in energy-rescaled Mino time
//! @param out state after the requested step; may alias y
//! @result Fills out; non-finite values must be rejected by the caller.
{
    // sample intermediate states to keep the fourth-order step self-contained.
    double k1[GEOD_N], k2[GEOD_N], k3[GEOD_N], k4[GEOD_N], tmp[GEOD_N];
    geodesic_priv_rhs(p, y, k1);
    for (int j=0; j<GEOD_N; ++j) tmp[j] = y[j]+h*k1[j]/2.0;
    geodesic_priv_rhs(p, tmp, k2);
    for (int j=0; j<GEOD_N; ++j) tmp[j] = y[j]+h*k2[j]/2.0;
    geodesic_priv_rhs(p, tmp, k3);
    for (int j=0; j<GEOD_N; ++j) tmp[j] = y[j]+h*k3[j];
    geodesic_priv_rhs(p, tmp, k4);
    // combine the stages into one update without modifying the input state.
    for (int j=0; j<GEOD_N; ++j)
        out[j] = y[j]+h*(k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])/6.0;
}

DEVICEFUNC static void geodesic_priv_two_steps(const geod_equations *p, const double y[GEOD_N], double h,
                      double out[GEOD_N])
//! Advance the geodesic with two equal RK4 half-steps.
//! - provides the refined solution used in the step-doubling error estimate
//! - also supplies consistent intermediate states during event localization
//! @param p immutable geodesic coefficients
//! @param y initial (u,v,m,w) state
//! @param h total Mino-time interval covered by both steps
//! @param out refined final state; may alias y
//! @result Fills out without modifying p.
{
    // two half-steps also provide an error estimate against one full step.
    double half[GEOD_N];
    geodesic_priv_rk4(p, y, h/2.0, half);
    geodesic_priv_rk4(p, half, h/2.0, out);
}

// locate a bracketed event by reintegrating from the start of the step.
DEVICEFUNC static double geodesic_priv_event_step(const geod_equations *p, const double y[GEOD_N], double h,
                         int component, double value)
//! Locate a bracketed coordinate event within one accepted integration step.
//! - bisects the interval and reintegrates from its initial state
//! - requires a sign change of y[component]-value over [0,h]
//! - the caller must ensure that this bracket contains only the desired event
//! @param p immutable geodesic coefficients
//! @param y state at the beginning of the bracket
//! @param h accepted step containing the event
//! @param component coordinate index, either inverse radius or cos(theta)
//! @param value coordinate value defining the surface
//! @result Mino-time offset in [0,h], refined through 48 bisections.
{
    // retain a sign-changing bracket while refining the event within this step.
    double lo=0.0, hi=h, trial[GEOD_N];
    for (int i=0; i<48; ++i) {
        double mid=(lo+hi)/2.0;
        geodesic_priv_two_steps(p, y, mid, trial);
        if ((trial[component]-value)*(y[component]-value) > 0.0) lo=mid;
        else hi=mid;
    }
    return (lo+hi)/2.0;
}

// bracket real roots of the paper's radial quartic between its stationary
// points. Direct quartic formulae lose the small roots when |E| is near
// zero and the large roots are many orders of magnitude farther away.
DEVICEFUNC static long double geodesic_priv_radial_potential(const geod_equations *p, long double r)
//! Evaluate the depressed radial quartic using extended-precision arithmetic.
//! - represents R(r)=r^4+A*r^2+2*B*r-C from Bursa (2017), Eq. (21)
//! - Horner evaluation reduces cancellation when bracketing small real roots
//! @param p coefficients of the energy-normalized radial potential
//! @param r Boyer-Lindquist radius in GM/c^2, including negative root candidates
//! @result Value of R(r); this routine does not test whether r is physically accessible.
{
    // evaluate the radial polynomial with extended precision near small roots.
    return ((r*r+p->A)*r+2.0L*p->B)*r-p->C;
}

DEVICEFUNC static int geodesic_priv_radial_roots(const geod_equations *p, double roots[4])
//! Find the real roots of the Kerr radial quartic without quartic-radical cancellation.
//! - uses the derivative cubic to split the real axis into monotonic intervals
//! - brackets and bisects each sign-changing interval in extended precision
//! - treats the exact radial Schwarzschild polynomial R=r^4 separately
//! @param p finite normalized polynomial coefficients
//! @param roots output real roots in ascending order; only returned entries are valid
//! @result Number of real roots (0,2,4), or -1 for a degenerate/failed calculation.
//! Exact repeated roots other than R=r^4 are deliberately reported as unsupported.
{
    // preserve the exact radial Schwarzschild limit despite its repeated root.
    if (p->A == 0.0 && p->B == 0.0 && p->C == 0.0) {
        for (int j=0; j<4; ++j) roots[j]=0.0;
        return 4;
    }
    // stationary points split the quartic into monotonic intervals.
    double cr[3], ci[3];
    int nc=cubic_eq(0.0,p->A/2.0,p->B/2.0,cr,ci);
    if (nc != 1 && nc != 3) return -1;
    long double knots[5];
    double bound=1.0+2.0*fmax(sqrt(fabs(p->A)),
                            fmax(cbrt(2.0*fabs(p->B)),pow(fabs(p->C),0.25)));
    knots[0]=-bound;
    // refine the stationary points before using them as root brackets.
    int n=1;
    for (int j=0; j<3; ++j) if (ci[j] == 0.0) {
        long double x=cr[j];
        for (int it=0; it<12; ++it) {
            long double f=(2.0L*x*x+p->A)*x+p->B;
            long double df=6.0L*x*x+p->A;
            if (df == 0.0L) break;
            x -= f/df;
        }
        if (!isfinite(x)) return -1;
        knots[n++]=x;
    }
    // sort the interior boundaries and enclose every possible real root.
    if (n != nc+1) return -1;
    for (int j=1; j<n; ++j) for (int k=j+1; k<n; ++k)
        if (knots[k] < knots[j]) {
            long double tmp=knots[j]; knots[j]=knots[k]; knots[k]=tmp;
        }
    knots[n++]=bound;
    // a sign change on a monotonic interval identifies exactly one real root.
    int nr=0;
    for (int j=0; j<n-1; ++j) {
        long double lo=knots[j], hi=knots[j+1];
        long double flo=geodesic_priv_radial_potential(p,lo), fhi=geodesic_priv_radial_potential(p,hi);
        // a precisely repeated root is a degenerate orbit, not a miss.
        if (flo == 0.0L || fhi == 0.0L) return -1;
        if ((flo > 0.0L) == (fhi > 0.0L)) continue;
        // bisect instead of subtracting nearly equal quartic radicals.
        for (int it=0; it<100; ++it) {
            long double mid=(lo+hi)/2.0L;
            long double f=geodesic_priv_radial_potential(p,mid);
            if ((f > 0.0L) == (flo > 0.0L)) lo=mid;
            else hi=mid;
        }
        roots[nr++]=(double)((lo+hi)/2.0L);
    }
    return nr;
}

// gL8 on smooth compact coordinates. CC has no real radial turning point.
DEVICEFUNC static double geodesic_priv_cc_integrand(geodesic *g, double x, int mode)
//! Smooth numerical integrands for the no-real-root (CC) radial branch.
//! - position and azimuth modes use u=1/r to include infinity at u=0
//! - the time mode uses x=log(r-r_h) to resolve horizon and distant endpoints
//! - the time mode contains only the radial part of Bursa (2017), Eq. (19)
//! @param g initialized CC geodesic
//! @param x inverse radius for position/azimuth, log(r-r_h) for time
//! @param mode 0 for position, 1 for radial coordinate time, 2 for radial azimuth
//! @result Integrand including the coordinate-change Jacobian; may be non-finite
//! if the caller supplies coordinates outside the exterior integration domain.
{
    // inverse radius keeps the integral to infinity on a finite interval.
    double a2=g->a*g->a, A=a2-g->l*g->l-g->q;
    double B=g->q+sqr(g->l-g->a), C=a2*g->q;
    if (mode!=1) {
        double value=1.0/sqrt(1.0+x*x*(A+x*(2.0*B-C*x)));
        // the radial azimuth rate is Eq. (20), expressed in inverse radius.
        if (mode==2) value *= (2.0*g->a*x-g->l*a2*x*x)/(1.0-2.0*x+a2*x*x);
        return value;
    }
    // logarithmic distance from the horizon resolves the time integrand at both ends.
    double rh=r_bh(g->a), dr=exp(x), r=rh+dr;
    double delta=dr*(r-(1.0-sqrt(1.0-a2)));
    double R=((r*r+A)*r+2.0*B)*r-C;
    return (r*r*(r*r+a2)+2.0*g->a*r*(g->a-g->l))/delta/sqrt(R)*dr;
}

DEVICEFUNC static double geodesic_priv_cc_integral(geodesic *g, double lo, double hi, int mode)
//! Integrate a CC radial position, time or azimuth integrand by Gauss-Legendre quadrature.
//! - doubles the number of eight-point panels until successive estimates agree
//! - uses a mixed absolute/relative tolerance of 1e-11 and at most 16384 panels
//! - preserves the sign of the integration bounds and uses no global state
//! @param g initialized CC geodesic
//! @param lo lower bound in inverse radius or log(r-r_h), according to mode
//! @param hi upper bound in the same coordinate
//! @param mode 0 for position, 1 for radial coordinate time, 2 for radial azimuth
//! @result Signed integral, zero for equal bounds, or NaN on numerical failure.
{
    // symmetric Gauss-Legendre nodes resolve each smooth subinterval.
    const double nodes[4]={0.1834346424956498,0.5255324099163290,
                           0.7966664774136267,0.9602898564975363};
    const double weights[4]={0.3626837833783620,0.3137066458778873,
                             0.2223810344533745,0.1012285362903763};
    // refine until consecutive integral estimates agree, or report non-convergence.
    if (lo == hi) return 0.0;
    double previous=NAN;
    for (int n=1; n<=16384; n*=2) {
        double h=(hi-lo)/n, sum=0.0;
        for (int j=0; j<n; ++j) {
            double mid=lo+(j+0.5)*h;
            for (int i=0; i<4; ++i)
                sum += h*0.5*weights[i]*(geodesic_priv_cc_integrand(g,mid-h*0.5*nodes[i],mode)
                                      +geodesic_priv_cc_integrand(g,mid+h*0.5*nodes[i],mode));
        }
        if (!isfinite(sum)) return NAN;
        if (n>1 && fabs(sum-previous)<1e-11*(1.0+fabs(sum))) return sum;
        previous=sum;
    }
    return NAN;
}

// q<0: m = sign(m) sqrt(m2p) dn(phase|mm). Fold at 2K and use
// the Jacobi epsilon function to integrate m^2 through any polar bounces.
DEVICEFUNC static double geodesic_priv_dn_phase(geodesic *g, double m)
//! Principal Jacobi dn phase for a vortical (q<0) polar coordinate.
//! - uses the positive root parameters already stored in the geodesic
//! - ignores the hemisphere sign; direction and hemisphere are tracked separately
//! - clamps turning-point roundoff before evaluating the inverse Jacobi function
//! @param g initialized vortical geodesic
//! @param m allowed cos(theta) coordinate
//! @result Dimensionless principal phase in [0,K(mm)].
{
    // clamp roundoff at polar turning points before inverting the dn phase.
    double z=(g->m2p-m*m)/(g->m2p+g->m2m);
    return jacobi_isn(sqrt(fmax(0.0,fmin(1.0,z))),g->mm);
}

DEVICEFUNC static double geodesic_priv_dn(geodesic *g, double phase)
//! Evaluate Jacobi dn after reducing an arbitrary phase to one period.
//! - permits negative phases and trajectories with repeated polar oscillations
//! - keeps the argument within the domain accepted by SIM5's Jacobi implementation
//! @param g initialized vortical geodesic supplying elliptic parameter mm
//! @param phase unwrapped dimensionless polar phase
//! @result dn(phase|mm); hemisphere and amplitude are applied by the caller.
{
    // fold into one period accepted by the Jacobi function implementation.
    double period=2.0*elliptic_k(g->mm);
    phase -= floor(phase/period)*period;
    return jacobi_dn(phase,g->mm);
}

DEVICEFUNC static double geodesic_priv_epsilon(geodesic *g, double phase)
//! Unwrapped Jacobi epsilon function for integrating dn-squared over polar motion.
//! - retains complete periods and reflects the remainder onto the principal branch
//! - uses the incomplete elliptic integral of the second kind for that remainder
//! @param g initialized vortical geodesic supplying elliptic parameter mm
//! @param phase arbitrary signed dimensionless phase
//! @result Integral of dn(u|mm)^2 from zero to phase, including complete periods.
{
    // retain complete periods so the polar time integral counts every oscillation.
    double K=elliptic_k(g->mm), E=elliptic_e(PI_half,g->mm);
    double n=floor(phase/(2.0*K)), x=phase-n*2.0*K;
    // reflect the remaining segment onto the principal incomplete-integral branch.
    double folded=fmin(x,2.0*K-x);
    double sn=jacobi_sn(folded,g->mm);
    double partial=elliptic_e(asin(fmax(-1.0,fmin(1.0,sn))),g->mm);
    return 2.0*n*E+(x>K ? 2.0*E-partial : partial);
}

DEVICEFUNC static double geodesic_priv_vortical_time(geodesic *g, double P1, double P2)
//! Polar contribution to coordinate travel time for a q<0 geodesic.
//! - integrates a^2*m^2 over the actual dn phase anchored by g->Tip
//! - uses absolute path parameters to avoid direction ambiguity at repeated endpoints
//! - does not include the radial contribution to the time delay
//! @param g vortical geodesic initialized from a source or from infinity
//! @param P1 positional integral at the first endpoint
//! @param P2 positional integral at the second endpoint
//! @result Nonnegative polar travel-time contribution in GM/c^3.
{
    // tip anchors the unwrapped dn phase at P=0; keep the actual phase,
    // since endpoints alone can be ambiguous over half-period intervals.
    double phase1=(g->Tip+P1)/g->mK, phase2=(g->Tip+P2)/g->mK;
    return sqr(g->a)*g->m2p*g->mK*fabs(geodesic_priv_epsilon(g,phase2)
                                    -geodesic_priv_epsilon(g,phase1));
}
//! \endcond

DEVICEFUNC
int geodesic_find_src_intersections(geodesic *g, double r, double m, double k[4],
    double target_r, geodesic_intersection *equator, geodesic_intersection *sphere, int *error)
//! Find the first future exterior equatorial and target-sphere intersections.
//! - uses the source geodesic initialized by geodesic_init_src() for the same r,m,k
//! - treats the equatorial plane as transparent and includes crossings after the sphere
//! - integrates separated motion with RK4 step doubling and locates events by bisection
//! - reports each surface independently; a captured photon may still cross the equator
//! @param g initialized geodesic; not modified by this routine
//! @param r source Boyer-Lindquist radius in GM/c^2, strictly outside the horizon
//! @param m source cos(theta), strictly between -1 and +1
//! @param k future-directed null momentum at the source, in Boyer-Lindquist components
//! @param target_r finite target sphere radius in GM/c^2, strictly greater than r
//! @param equator output first crossing; found=0 and NaN coordinates mean no crossing
//! @param sphere output first sphere arrival; found=0 and NaN coordinates mean no arrival
//! @param error optional output: GD_OK, GD_ERROR_ARGUMENT, or GD_ERROR_NUMERICAL
//! @result TRUE when both intersection questions are resolved, FALSE on failure.
//! Ignore both outputs on FALSE. Returned P uses geodesic_P_int()'s branch convention;
//! supported endpoint pairs can be passed to geodesic_timedelay(). The local tolerance
//! is 1e-11 and the step-attempt limit is 200000; neither guarantees relative radius
//! accuracy for nearly degenerate or extremely distant crossings. No global state is used.
{
    // distinguish invalid input from a physical miss and initialize both outputs.
    if (error) *error=GD_ERROR_ARGUMENT;
    if (!g || !k || !equator || !sphere) return FALSE;
    equator->found=sphere->found=0;
    equator->P=equator->r=equator->m=sphere->P=sphere->r=sphere->m=NAN;
    if (!isfinite(r) || !isfinite(m) || !isfinite(target_r) ||
        r<=r_bh(g->a) || fabs(m)>=1.0 || target_r<=r) return FALSE;
    // derive the separated equations from the already initialized geodesic.
    if (error) *error=GD_ERROR_NUMERICAL;
    geod_equations p={g->a,g->l,g->q,sqr(g->a)-sqr(g->l)-g->q,
                      g->q+sqr(g->l-g->a),sqr(g->a)*g->q};
    sim5metric metric;
    kerr_metric(g->a,r,m,&metric);
    double E=-(metric.g00*k[0]+metric.g03*k[3]), energy=fabs(E);
    if (!isfinite(energy) || energy==0.0) return FALSE;
    // exterior turning points determine whether the target sphere is reachable.
    double horizon=r_bh(g->a), roots[4]={creal(g->r1),creal(g->r2),creal(g->r3),creal(g->r4)};
    int inner=0,outer=0;
    for (int j=0; j<g->nrr; ++j) {
        if (roots[j]>horizon && roots[j]<r) inner=1;
        if (roots[j]>r && roots[j]<target_r) outer=1;
    }
    int reaches=E>0.0 && !outer && (k[1]>=0.0 || inner);
    int cross_known=g->q<=0.0, target_known=!reaches;
    // preserve the library parameter convention on each radial branch.
    int radial=(g->a==0.0 && g->l==0.0 && g->q==0.0);
    int branch=(g->type==GEOD_TYPE_RR_BH) ? k[1]<0.0 : k[1]>=0.0;
    double P0=geodesic_P_int(g,r,branch), direction=1.0;
    if (g->type==GEOD_TYPE_RC || g->type==GEOD_TYPE_CC || radial)
        direction=k[1]>=0.0 ? -1.0 : 1.0;
    if (!isfinite(P0)) return FALSE;
    // absolute energy keeps negative-energy rays future-directed inside the ergosphere.
    double sigma=r*r+g->a*g->a*m*m;
    double y[GEOD_N]={1.0/r,-sigma*k[1]/(energy*r*r),m,
                     -sqrt(1.0-m*m)*sigma*k[2]/energy};
    double elapsed=0.0;
    // a source on the equator already supplies its first intersection.
    if (m==0.0) {
        equator->found=1; equator->P=P0; equator->r=r; equator->m=0.0;
        cross_known=1;
    }
    // limit phase advance so a single step cannot skip a polar oscillation.
    double hmax=0.05/sqrt(1.0+fabs(p.A)+fabs(p.B)+fabs(p.C)), h=hmax;
    for (int step=0; step<200000; ++step) {
        // stop when both questions are answered, retaining any genuine misses.
        if (cross_known && target_known) { if (error) *error=GD_OK; return TRUE; }
        if (h<1e-15) return FALSE;
        // reject inaccurate steps before accepting positions or testing intersections.
        double full[GEOD_N],fine[GEOD_N],err=0.0;
        geodesic_priv_rk4(&p,y,h,full);
        geodesic_priv_two_steps(&p,y,h,fine);
        for (int j=0; j<GEOD_N; ++j) {
            if (!isfinite(full[j]) || !isfinite(fine[j])) return FALSE;
            err=fmax(err,fabs(fine[j]-full[j])/(15.0*1e-11*(1.0+fabs(fine[j]))));
        }
        if (err>1.0) { h*=0.5; continue; }
        // process the earliest event so an interior crossing cannot precede capture.
        int event=-1;
        double event_h=h,eh;
        if (fine[GEOD_U]>=1.0/horizon) {
            event_h=geodesic_priv_event_step(&p,y,h,GEOD_U,1.0/horizon); event=0;
        }
        if (fine[GEOD_U]<=0.0) {
            eh=geodesic_priv_event_step(&p,y,h,GEOD_U,0.0);
            if (eh<=event_h) { event_h=eh; event=1; }
        }
        if (!cross_known && fine[GEOD_M]*m<=0.0) {
            eh=geodesic_priv_event_step(&p,y,h,GEOD_M,0.0);
            if (eh<event_h) { event_h=eh; event=2; }
        }
        if (!target_known && fine[GEOD_U]<=1.0/target_r) {
            eh=geodesic_priv_event_step(&p,y,h,GEOD_U,1.0/target_r);
            if (eh<event_h) { event_h=eh; event=3; }
        }
        // advance only to the selected event, not to the discarded step endpoint.
        if (event>=0) geodesic_priv_two_steps(&p,y,event_h,fine);
        elapsed+=event_h;
        for (int j=0; j<GEOD_N; ++j) y[j]=fine[j];
        // check the predicted radial fate before declaring the remaining surface a miss.
        if (event==0 || event==1) {
            if ((event==0 && reaches) || (event==1 && !sphere->found)) return FALSE;
            if (error) *error=GD_OK;
            return TRUE;
        }
        // return coordinates together with the parameter needed by the time-delay API.
        if (event==2 || event==3) {
            geodesic_intersection *hit=event==2 ? equator : sphere;
            hit->found=1; hit->P=P0+direction*elapsed;
            hit->r=event==2 ? 1.0/y[GEOD_U] : target_r;
            hit->m=event==2 ? 0.0 : y[GEOD_M];
            if (event==2) cross_known=1;
            else target_known=1;
        }
        // recover larger steps only when the error estimate leaves enough margin.
        if (err<0.03) h=fmin(2.0*h,hmax);
    }
    return FALSE;
}
