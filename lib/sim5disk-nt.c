//************************************************************************
//    sim5disk-nt.c
//------------------------------------------------------------------------
//    Date   : 2.10.2014
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2014 Michal Bursa
//************************************************************************





static float bh_mass    = 10.0;
static float bh_spin    = 0.0;
static float disk_mdot  = 0.1;
static float disk_rms   = 6.0;
static float disk_alpha = 0.1;
static int   options    = 0;

int disk_nt_setup(double M, double a, double mdot_or_L, double alpha, int _options)
{
    bh_mass    = M;
    bh_spin    = a;
    disk_rms   = disk_nt_r_min();
    disk_alpha = alpha;
    options    = _options;
    if (options & DISK_NT_OPTION_LUMINOSITY) {
        disk_mdot = mdot_or_L;
        double find_mdot_for_luminosity(double L0);
        disk_mdot = find_mdot_for_luminosity(mdot_or_L);
        fprintf(stderr,"(disk-nt) final mdot for L=%.5f: mdot=%.5f (%.6e 10^18 g/s)\n", mdot_or_L, disk_mdot, disk_mdot*bh_mass*Mdot_Edd/1e18);
    }
    else {
        disk_mdot = mdot_or_L;
        fprintf(stderr,"(disk-nt) mdot set to %.5f: (%.6e 10^18 g/s; lum=%.5f)\n", disk_mdot, disk_mdot*bh_mass*Mdot_Edd/1e18, disk_nt_lum());
    }
    return 0;
}



void disk_nt_finish()
{
}



double disk_nt_r_min()
{
    double a = bh_spin;
    double z1,z2,r0;
    double sga = (a>=0.0) ? +1. : -1.;
    z1 = 1.+pow(1.-a*a, 1./3.)*(pow(1.+a, 1./3.)+pow(1.-a, 1./3.));
    z2 = sqrt(3.*a*a+z1*z1);
    r0 = 3.+z2-sga*sqrt((3.-z1)*(3.+z1+2.*z2));
    return r0+1e-3;
}



double disk_nt_flux(double r)
{
    if (r <= disk_rms) return 1e-40;
    double a = bh_spin;
    double x=sqrt(r);
    double x0,x1,x2,x3;
    x0=sqrt(disk_rms);
    x1=+2.*cos(1./3.*acos(a)-M_PI/3.);
    x2=+2.*cos(1./3.*acos(a)+M_PI/3.);
    x3=-2.*cos(1./3.*acos(a));
    double f0,f1,f2,f3,F;
    f0=x-x0-1.5*a*log(x/x0);
    f1=3.*sqr(x1-a)/(x1*(x1-x2)*(x1-x3))*log((x-x1)/(x0-x1));
    f2=3.*sqr(x2-a)/(x2*(x2-x1)*(x2-x3))*log((x-x2)/(x0-x2));
    f3=3.*sqr(x3-a)/(x3*(x3-x1)*(x3-x2))*log((x-x3)/(x0-x3));
    F = 1./(4.*M_PI*r) * 1.5/(x*x*(x*x*x-3.*x+2.*a)) * (f0-f1-f2-f3);
    return 9.1721376255e+28 * F * disk_mdot/bh_mass;
}



double disk_nt_lum()
{
    const float disk_rmax = 100000.;
    double r1;
    int nflux;

    nflux=0;
    for (r1=disk_rms; r1<disk_rmax; r1*=1.005) nflux++;

    double* flux[2];
    flux[0] = (double*) calloc(2*nflux, sizeof(double));
    flux[1] = (double*) calloc(2*nflux, sizeof(double));
    nflux=0;

    for(r1=disk_rms;r1<disk_rmax;r1*=1.005) {
        flux[0][nflux] = r1*bh_mass*grav_radius; // radius in [cm]
        flux[1][nflux] = 2.0*disk_nt_flux(r1); // flux in [erg cm-2 s-1];
        nflux++;
    }

    gsl_spline *splineF;
    gsl_interp_accel *accF;
    splineF = gsl_spline_alloc(gsl_interp_linear, nflux);
    accF    = gsl_interp_accel_alloc();
    gsl_spline_init(splineF, flux[0], flux[1], nflux);



    double func_luminosity(double x, void* params)
    {
        return 2.*M_PI*x*gsl_spline_eval(splineF,x,accF);
    }

    gsl_integration_workspace *iws = gsl_integration_workspace_alloc(1000);
    gsl_function F;
    F.function = &func_luminosity;
    F.params = NULL;
    double L, err;
    gsl_integration_qag(&F, flux[0][0], flux[0][nflux-1], 0, 1e-3, 1000, GSL_INTEG_GAUSS15, iws, &L, &err);
    gsl_integration_workspace_free(iws);
    free(flux[0]);
    free(flux[1]);
    return L/(L_Edd*bh_mass);
}



double disk_nt_mdot()
{
    return disk_mdot;
}



double disk_nt_temp(double r)
{
    return sqrt4(disk_nt_flux(r)/sb_sigma);
}



double disk_nt_sigma(double r)
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



double disk_nt_ell(double r)
{
    double a = bh_spin;
    r = max(disk_rms, r);
    return (r*r-2.*a*sqrt(r)+a*a) / (sqrt(r)*r-2.*sqrt(r)+a);
}



double disk_nt_vr(double r)
{
    return 0.0;
}



double disk_nt_h(double r)
{
    return 0.0;
}



double disk_nt_dhdr(double r)
{
    return 0.0;
}



void disk_nt_dump()
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



double find_mdot_for_luminosity(double L0) {
    double L;

    double fce(double xmdot) {
        disk_mdot = xmdot;
        return L0-disk_nt_lum();
    }

    int res = rtbis(0.0, 100.0, 1e-6, fce, &L);
    return (res) ? L : 0.0;
}   


