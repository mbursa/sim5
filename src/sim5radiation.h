//************************************************************************
//    sim5radiation.h
//------------------------------------------------------------------------
//    Date   : 2.10.2014
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2014 Michal Bursa
//************************************************************************


#ifndef _SIM5RADIATION_H
#define _SIM5RADIATION_H


typedef struct stokes_params {
  double i;
  double q;
  double u;
  double v;
} stokes_params;


static const stokes_params stokes_null = {0.0, 0.0, 0.0, 0.0};



DEVICEFUNC double blackbody_Iv(double T, double hardf, double cos_mu, double E);
DEVICEFUNC void blackbody(double T, double hardf, double cos_mu, double E[], double Iv[], int en_bins);
DEVICEFUNC INLINE double blackbody_photons(double T, double hardf, double cos_mu, double E);
DEVICEFUNC double blackbody_photons_total(double T, double hardf, double cos_mu);
DEVICEFUNC double blackbody_photon_energy_random(double T);


#endif
