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


void blackbody(double T, double hardf, double cos_mu, double E[], double Iv[], int en_bins);
double blackbody_Iv(double T, double hardf, double cos_mu, double E);



#endif
