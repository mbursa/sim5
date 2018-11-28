//************************************************************************
//    SIM5 library
//    sim5radiation.c - ddd
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute of the Czech Academy of Sciences
//************************************************************************


//! \file sim5radiation.c
//! Radiative processes routines
//! 
//! Provides routines for radiative processes.






//-----------------------------------------------------------------------------------
// blackbody
//-----------------------------------------------------------------------------------


DEVICEFUNC
double blackbody_Iv(double T, double hardf, double cos_mu, double E)
//! Specific radiance of black-body radiation.
//!
//! Gives specific intensity \f\frac{dE}{dt dA d\Omega dE}\f] at energy `E` of radiation
//! of a black body of temperatute `T`. Assumes either isotropic or limb-darkened emission and
//! a color correction of effective temperature by factor `hardf`.
//!
//! @param T temperature [K]
//! @param hardf hardening factor (effective temperature correction)
//! @param cos_mu cosine of emission direction with respect to the normal to the emission surface;
//!               set cos_mu>=0 for limb-darkened emission and cos_mu<0 for isotropic emission
//! @param E energy [keV]
//!
//! @result specific intensity in units [erg cm^-2 s^-1 keV^-1 srad^-1]
{
    if (T<=0.0) return 0.0;
    // make spectrum using black-body formula with hardening factor & limb darkening
    double limbf = (cos_mu>=0.0) ? 0.5+0.75*cos_mu : 1.0;
    double freq  = kev2freq*E;
    return limbf * 2.0*planck_h*sqr3(freq)/sqr(speed_of_light)/sqr4(hardf) / expm1((planck_h*freq)/(boltzmann_k*hardf*T)) * (1./freq2kev);
    // freq2kev factor converts frequency to keV in denominator
}



DEVICEFUNC
void blackbody(double T, double hardf, double cos_mu, double E[], double Iv[], int N)
//! Specific radiance of black-body radiation.
//!
//! Gives specific intensity \f$\frac{dE}{dt dA d\Omega dE}\f$ at energy `E` of radiation
//! of a black body of temperatute `T`. Assumes either isotropic or limb-darkened emission and
//! a color correction of effective temperature by factor `hardf`.
//!
//! @param T temperature [K]
//! @param hardf hardening factor (effective temperature correction)
//! @param cos_mu cosine of emission direction with respect to the normal to the emission surface;
//!               set cos_mu>=0 for limb-darkened emission and cos_mu<0 for isotropic emission
//! @param E array of energies (input) [keV]
//! @param Iv array of specific intensities (output) [erg cm^-2 s^-1 keV^-1 srad^-1]
//! @param N dimension of the E[] and Iv[] arrays
//!
//! @result `Iv[]` array contains specific intensities for given energies
{
    if (T<=0.0) return;

    // make spectrum using black-body formula with hardening factor & limb darkening
    int i;
    double limbf = (cos_mu>=0.0) ? 0.5+0.75*cos_mu : 1.0;
    double BB1 = limbf * 2.0*planck_h/sqr(speed_of_light)/sqr4(hardf)*sqr4(kev2freq);
    double BB2 = (planck_h*kev2freq)/(boltzmann_k*hardf*T);
    for (i=0; i<N; i++) Iv[i] = BB1*sqr3(E[i])/expm1(BB2*E[i]);
}



DEVICEFUNC INLINE
double blackbody_photons(double T, double hardf, double cos_mu, double E)
//! Specific photon intensity of black-body radiation.
//!
//! Same as `blackbody_Iv()` function, except it gives specific photon intensity \f$\frac{dN}{dt dA d\Omega dE}\f$.
//!
//! @result specific intensity in units [photons cm^-2 s^-1 keV^-1 srad^-1]
{
    return blackbody_Iv(T, hardf, cos_mu, E) / (E*kev2erg);
}



DEVICEFUNC
double blackbody_photons_total(double T, double hardf)
//! Number of photons that are emitted by a black-body surface.
//!
//! Gives total number of photons \f$\frac{dN}{dt dA d\Omega}\f$ of all energies that are
//! emitted form a unit surface of a black body of temperatute `T` into a unit solid angle.
//! Assumes either isotropic or limb-darkened emission and a color correction of
//! effective temperature by factor `hardf`.
//!
//! @param T temperature [K]
//! @param hardf hardening factor (effective temperature correction)
//! @param cos_mu cosine of emission direction with respect to the normal to the emission surface;
//!               set cos_mu>=0 for limb-darkened emission and cos_mu<0 for isotropic emission
//!
//! @result number of photons [photons cm^-2 s^-1 srad^-1]
{
    // 4.808227612 = 4*Zeta(3)
    return M_PI * 4.808227612 * sqr3(T) * sqr3(boltzmann_k) / sqr3(planck_h) / speed_of_light2 / hardf;
}



DEVICEFUNC
double blackbody_photon_energy_random(double T)
//! Draws a random photon energy that follows Planck distribution.
//!
//! Picks a random photon energy of black-body radiation according to
//! Planck energy distribution at given temperature.
//! For derivation see http://arxiv.org/abs/1307.3635, part 3.3.1.
//!
//! @param T temperature of the distribution [K] (including hardening factor)
//!
//! @result photon energy in [keV]
{
    double u1 = sim5urand();
    double u2 = sim5urand();
    double u3 = sim5urand();
    double u4 = sim5urand();
    int m;
    for (m=1; ; m++) {
        double sum_j = 0.0;
        for (int i=1; i<=m; i++) sum_j += 1.0/(double)(i*i*i);
        if (1.202*u1 < sum_j) break;
    }
    return boltzmann_k*T * (-log(u2*u3*u4)) / (double)(m) * erg2kev;
}



