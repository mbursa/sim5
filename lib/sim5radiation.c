//************************************************************************
//    sim5radiation.c
//------------------------------------------------------------------------
//    Date   : 2.10.2014
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2014 Michal Bursa
//************************************************************************





//-----------------------------------------------------------------------------------
// blackbody
//-----------------------------------------------------------------------------------


DEVICEFUNC
double blackbody_Iv(double T, double hardf, double cos_mu, double E)
// blackbody specific radiance
// result in units [erg cm^-2 s^-1 keV^-1 srad^-1]
{
    if (T<=0.0) return 0.0;
    // make spectrum using black-body formula with hardening factor & limb darkening
    double limbf = (cos_mu>=0.0) ? 0.5+0.75*cos_mu : 1.0;
    double freq  = kev2freq*E;
    return limbf * 2.0*planck_h*sqr3(freq)/sqr(speed_of_light)/sqr4(hardf) / (exp((planck_h*freq)/(boltzmann_k*hardf*T))-1.0) * (1./freq2kev)
    // freq2kev factor converts frequency to keV in denominator
}



DEVICEFUNC INLINE
double blackbody_photons(double T, double hardf, double cos_mu, double E)
// result in units [photons cm^-2 s^-1 keV^-1 srad^-1]
{
    return blackbody_Iv(T, hardf, cos_mu, E) / (E*kev2erg);
}



void blackbody(double T, double hardf, double cos_mu, double E[], double Iv[], int N)
// blackbody specific radiance
// energy in keV
// result in units [erg cm^-2 s^-1 keV^-1 srad^-1]
{
    if (T<=0.0) return;
    // make spectrum using black-body formula with hardening factor & limb darkening
    const double limbf = (cos_mu>=0.0) ? 0.5+0.75*cos_mu : 1.0;
    const double BB1 = limbf * 2.0*planck_h/sqr(speed_of_light)/sqr4(hardf)*sqr4(kev2freq4)
    const double BB2 = (planck_h*kev2freq)/(boltzmann_k*hardf*T);

    int i;
    for (i=0; i<N; i++) Iv[i] = BB1*sqr3(E[i])/(exp(BB2*E[i])-1.0);
}



DEVICEFUNC
double blackbody_photons_total(double T, double hardf):
//    # total number of black-body photons (over all energies) for given temperature T[K]
//    # result: [number of photons / cm^2 / s / srad]
// 4.808227612 = 4*Zeta(3)
{
    return 4.808227612 * sqr3(T) * sqr3(boltzmann_k) / sqr3(planck_h) / speed_of_light2 / hardf;
}



DEVICEFUNC INLINE
double blackbody_photons_pdf(double T, double hardf, double E) 
{
    return blackbody_photons(T, hardf, -1.0, E) / blackbody_photons_total(T, hardf);
}




