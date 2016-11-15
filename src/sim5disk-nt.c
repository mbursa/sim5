//************************************************************************
//    sim5disk-nt.c
//------------------------------------------------------------------------
//    Date   : 2.10.2014
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2014 Michal Bursa
//************************************************************************

//! \file sim5disk-nt.c
//! Thin disk routines
//! 
//! Provides routines for Novikov-Thorne thin disk model



// TODO: these should be made thread local
static float bh_mass    = 10.0;
static float bh_spin    = 0.0;
static float disk_mdot  = 0.1;
static float disk_rms   = 6.0;
static float disk_alpha = 0.1;
static int   options    = 0;



DEVICEFUNC
int disk_nt_setup(double M, double a, double mdot_or_L, double alpha, int _options)
//! Set up a relativistic (Novikov-Thorne) model of thin disk.
//! The disk can be set up using either mass accretion rate or luminosity. In case of Mdot, this is 
//! specified as a ratio of actual accretion rate in grams per second to the Eddington mass acretion rate 
//! corresponding to the given mass. The Eddington mass accretion rate constant is declared in sim5const.h.
//! In case of luminosity, the model sets accretion rate such that the integrated disk luminosity matches 
//! the given value in ergs/sec relative to the Eddington luminosity again declared in sim5const.h.
//!
//! @param M mass of central BH [M_sun]
//! @param a spin of central BH [0..1]
//! @param mdot_or_L mass accretion rate (default) or luminosity (both in eddington units; see sim5const.h)
//! @param alpha viscosity parameter
//! @param options optional switches (default=0; switches can be combined with `+` operator)
//!        * DISK_NT_OPTION_LUMINOSITY: `mdot_or_L` parameter is interpreted as luminosity
//!
//! @result A status code (currently zero)
{
    bh_mass    = M;
    bh_spin    = a;
    disk_rms   = disk_nt_r_min();
    disk_alpha = alpha;
    options    = _options;
    if (options & DISK_NT_OPTION_LUMINOSITY) {
        disk_mdot = mdot_or_L;
        double disk_nt_find_mdot_for_luminosity(double L0);
        disk_mdot = disk_nt_find_mdot_for_luminosity(mdot_or_L);
        //fprintf(stderr,"(disk-nt) final mdot for L=%.5f: mdot=%.5f (%.6e 10^18 g/s)\n", mdot_or_L, disk_mdot, disk_mdot*bh_mass*Mdot_Edd/1e18);
    }
    else {
        disk_mdot = mdot_or_L;
        //fprintf(stderr,"(disk-nt) mdot set to %.5f: (%.6e 10^18 g/s; lum=%.5f)\n", disk_mdot, disk_mdot*bh_mass*Mdot_Edd/1e18, disk_nt_lum());
    }
    return 0;
}



DEVICEFUNC
void disk_nt_finish()
//! Finalize the disk model.
//! Cleans up and frees memory.
{
}



DEVICEFUNC
double disk_nt_r_min()
//! Minimal radius of the disk (disk inner edge).
//! Provides minimal value for radius for which the functions provide valid results. 
//! For NT disk, this corresponds to the radius of marginally stable orbit (r_ms), where there is zero torque in the fluid.
//!
//! @result Radius of disk inner edge [GM/c2]
{
    double a = bh_spin;
    double z1,z2,r0;
    double sga = (a>=0.0) ? +1. : -1.;
    z1 = 1.+pow(1.-a*a, 1./3.)*(pow(1.+a, 1./3.)+pow(1.-a, 1./3.));
    z2 = sqrt(3.*a*a+z1*z1);
    r0 = 3.+z2-sga*sqrt((3.-z1)*(3.+z1+2.*z2));
    return r0+1e-3;
}



DEVICEFUNC
double disk_nt_flux(double r)
//! Local flux from one side of the disk.
//! Provides radial radiation flux dependence for Novikov-Thorne accretion disk.
//! Formulae based on Page&Thorne(1974) http://adsabs.harvard.edu/abs/1974ApJ...191..499P
//!
//! Note the retuned flux is local flux, i.e. flux measured by observer that is at rest with respect to the fluid.
//!
//! @param r radius of emission [GM/c2]
//!
//! @result Total outgoing flux from unit area on one side of the disk [erg cm-2 s-1]. 
{
    if (r <= disk_rms) return 0.0;
    double a = bh_spin;
    double x=sqrt(r);
    double x0,x1,x2,x3;
    x0=sqrt(disk_rms);
    x1=+2.*cos(1./3.*acos(a)-M_PI/3.);
    x2=+2.*cos(1./3.*acos(a)+M_PI/3.);
    x3=-2.*cos(1./3.*acos(a));
    double f0,f1,f2,f3,F;
    // PT74 (eq.15n)
    f0=x-x0-1.5*a*log(x/x0);
    f1=3.*sqr(x1-a)/(x1*(x1-x2)*(x1-x3))*log((x-x1)/(x0-x1));
    f2=3.*sqr(x2-a)/(x2*(x2-x1)*(x2-x3))*log((x-x2)/(x0-x2));
    f3=3.*sqr(x3-a)/(x3*(x3-x1)*(x3-x2))*log((x-x3)/(x0-x3));
    F = 1./(4.*M_PI*r) * 1.5/(x*x*(x*x*x-3.*x+2.*a)) * (f0-f1-f2-f3);

    // in Newtonian limit flux is F = 3/(8pi) * G*M*Mdot/r^3 = 3*x^3/(8pi) G^-2*Mdot*M^-2*c^6
    // where x=r/(GM/c2), m=M/M_sun, mdot=Mdot/(M/M_sun*Mdot_Edd), r=GM/c2*x; Mdot_Edd is Eddington acc rate for BH of
    // F ~ mdot/m * Mdot_Edd * G^-2 * c^6 * M_sun^-2 [kg s^-3] = mdot/m * 9.172138e+28 [erg cm^-2 s^-1]
    // (Mdot_Edd*gram2kg)/sqr(grav_const)*sqr3(speed_of_light2)/sqr(solar_mass)*joule2erg/msq2cmsq = 9.1721376255e+28

    // the result shall be adjusted for actuall BH mass and acc rate following the scaling
    // F~mdot/m, where m=M/M_sun and mdot=Mdot/(M/M_sun*Mdot_Edd)

    return 9.1721376255e+28 * F * disk_mdot/bh_mass; // [erg cm^-2 s^-1]
}



DEVICEFUNC
double disk_nt_lum()
//! Gets total disk luminosity.
//! Luminosity is obtained by integrating local flux over the surface area of the disk 
//! going into the whole sky (4pi solid angle). 
//! The integration makes a proper transformation of local flux to coordinate frame, but 
//! ignores other relativistic effects like light bending.
//!
//! \f[L = 2 * 2\pi \int F(r) r dr\f]
//!
//! @result Total disk luminosity of both surfaces [erg s-1]
{
    const float disk_rmax = 1e5;

    // integrate disk luminosity from r_ms to disk_rmax rg
    // - the integration uses 'logarithmic rule': L = \int f(x) dx \int f(x)*x d(log(x))

    double func_luminosity(double log_r)
    {
        double r = exp(log_r);
        // calculate U_t
        double gtt = -1. + 2./r;
        double gtf = -2.*bh_spin/r;
        double gff = sqr(r) + sqr(bh_spin) + 2.*sqr(bh_spin)/r;
        double Omega = 1./(bh_spin + pow(r,1.5));
        double U_t = sqrt(-1.0/(gtt + 2.*Omega*gtf + sqr(Omega)*gff)) * (gtt + Omega*gtf);
        double F = disk_nt_flux(r);
        // dL = 2pi*r*F(r) dr, extra r comes from log integration
        return 2.*M_PI*r*2.0*(-U_t)*F * r;
    }

    double L = integrate_simpson(func_luminosity, log(disk_rms), log(disk_rmax), 1e-5);

    // fix units to erg/s
    L *= sqr(bh_mass*grav_radius);

    return L/(L_Edd*bh_mass);
}



DEVICEFUNC
double disk_nt_mdot()
//! Mass accretion rate.
//! Returns mass accretion rate in Eddington units. See `disk_nt_setup()` for details.
//!
//! @result Mass accretion rate in Eddington units.
{
    return disk_mdot;
}



DEVICEFUNC
double disk_nt_temp(double r)
//! Effective temperature.
//! Returns disk effective temperature at given radius. The temperature value corresponds to 
//! total outgoing flux at thar radius through the Steffan-Boltzmann law.
//!
//! @param r radius (measured in equatorial plane) [rg]
//!
//! @result Effective surface temperatute in [K].
{
    return sqrt4(disk_nt_flux(r)/sb_sigma);
}



DEVICEFUNC
double disk_nt_sigma(double r)
//! Column density.
//! Returns midplane column density of the fluid at given radius for the first two zones 
//! according to formulae from Chandrasekhar book.
//!
//! \f[ \Sigma = \int_0^H \rho dz \f]
//!
//! @param r radius (measured in equatorial plane) [rg]
//!
//! @result Midplane column density in [g/cm2].
{
    if (r < disk_rms) return 0.0;
    double a = bh_spin;

    double x=sqrt(r);
    double x0,x1,x2,x3;
    x0=sqrt(disk_rms);
    x1=+2.*cos(1./3.*acos(a)-M_PI/3.);
    x2=+2.*cos(1./3.*acos(a)+M_PI/3.);
    x3=-2.*cos(1./3.*acos(a));

    double xA, xB, xC, xD, xE, xL;
    xA = 1. + sqr(a)/sqr(r) + 2.*sqr(a)/sqr3(r);
    xB = 1. + a/(x*x*x);
    xC = 1. - 3./(x*x) + 2.*a/(x*x*x);
    xD = 1. - 2./r + sqr(a)/sqr(r);
    xE = 1. + 4.*sqr(a)/sqr(r) - 4.*sqr(a)/sqr3(r) + 3.*sqr4(a)/sqr4(r);

    double f0, f1, f2, f3;
    f0=x-x0-1.5*a*log(x/x0);
    f1=3.*(x1-a)*(x1-a)/(x1*(x1-x2)*(x1-x3))*log((x-x1)/(x0-x1));
    f2=3.*(x2-a)*(x2-a)/(x2*(x2-x1)*(x2-x3))*log((x-x2)/(x0-x2));
    f3=3.*(x3-a)*(x3-a)/(x3*(x3-x2)*(x3-x1))*log((x-x3)/(x0-x3));
    xL = (1.+a/(x*x*x))/sqrt(1.-3./(x*x)+2.*a/(x*x*x))/x * (f0-f1-f2-f3);

    double xMdot = disk_mdot*bh_mass*Mdot_Edd/1e17;
    double r_im = 40.*(pow(disk_alpha,2./21.)/pow(bh_mass/3.,2./3.)*pow(xMdot,16./20.)) * pow(xA,20./21.) *
    pow(xB,-36./21.) * pow(xD,-8./21.) * pow(xE,-10./21.) * pow(xL,16./21.);

    double Sigma;
    if (r < r_im)
        Sigma = 20. * (bh_mass/3.)/xMdot/disk_alpha * sqrt(r*r*r) * 1./(xA*xA) * pow(xB,3.) * sqrt(xC) * xE * 1./xL;
    else {
        Sigma = 5e4 * pow(bh_mass/3.,-2./5.)*pow(xMdot,3./5.)*pow(disk_alpha,-4./5.) * pow(r,-3./5.) * pow(xB,-4./5.) * sqrt(xC) * pow(xD,-4./5.) * pow(xL,3./5.);
    }

    return Sigma;
}



DEVICEFUNC
double disk_nt_ell(double r)
//! Specific angular momentum.
//! Returns specific angular momentum of the fluid at given radius. 
//!
//! @param r radius (measured in equatorial plane) [rg]
//!
//! @result Specific angular momentum in [g.u.].
{
    double a = bh_spin;
    r = max(disk_rms, r);
    return (r*r-2.*a*sqrt(r)+a*a) / (sqrt(r)*r-2.*sqrt(r)+a);
}



DEVICEFUNC
double disk_nt_vr(double r)
//! Radial velocity.
//! Returns bulk radial velocity of the fluid at given radius. 
//! (For thin disk approximation, this is always zero.)
//!
//! @param r radius (measured in equatorial plane) [rg]
//!
//! @result Radial velocity in [speed_of_light].
{
    return 0.0;
}



DEVICEFUNC
double disk_nt_h(double r)
//! Surface height.
//! Returns height of the surface of the disk (effective photosphere) above midplane at given radius. 
//! (For thin disk approximation, this is always zero.)
//!
//! @param r radius (measured in equatorial plane) [rg]
//!
//! @result Radial velocity in [speed_of_light].
{
    return 0.0;
}



DEVICEFUNC
double disk_nt_dhdr(double r)
//! Derivative of surface height.
//! Returns surface profile as derivative \f$dH/dR\f$ of its height above midplane at given radius. 
//! (For thin disk approximation, this is always zero.)
//!
//! @param r radius (measured in equatorial plane) [rg]
//!
//! @result Derivative of surface height.
{
    return 0.0;
}



DEVICEFUNC
void disk_nt_dump()
//! Dump function printing disk provile.
//! The function prints to stdout profile of all quantities as a function fo radius from r_ms to some outer radius (~2000 rg).
{
    const float disk_rmax = 2000.;
    printf("# (sim5disk-nt) dump\n");
    printf("#-------------------------------------------\n");
    printf("# M        = %.4f\n", bh_mass);
    printf("# a        = %.4f\n", bh_spin);
    printf("# rmin     = %.4f\n", disk_rms);
    printf("# rmax     = %.4f\n", disk_rmax);
    printf("# alpha    = %.4f\n", disk_alpha);
    printf("# options  = %d\n",   options);
    printf("# L        = %e\n", disk_nt_lum());
    printf("# mdot     = %e\n", disk_nt_mdot());
    printf("#-------------------------------------------\n");
    printf("# r   flux   sigma   ell   vr   H   dH/dr\n");
    printf("#-------------------------------------------\n");

    double r;
    for (r=disk_rms; r<disk_rmax; r*=1.05) {
        printf(
            "%e  %e  %e  %e  %e  %e  %e\n",
            r,
            disk_nt_flux(r),
            disk_nt_sigma(r),
            disk_nt_ell(r),
            disk_nt_vr(r),
            disk_nt_h(r),
            disk_nt_dhdr(r)
        );
    }
}



// private routine to iteratively find mdot that corresponds to given luminosity
DEVICEFUNC
double disk_nt_find_mdot_for_luminosity(double L0) {
    double L;

    double fce(double xmdot) {
        disk_mdot = xmdot;
        return L0-disk_nt_lum();
    }

    int res = rtbis(0.0, 100.0, 1e-6, fce, &L);
    return (res) ? L : 0.0;
}


