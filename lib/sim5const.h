//************************************************************************
//    SIM5 library
//    sim5const.h - fundamental physical constants and unit conversions
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************


#ifndef _SIM5CONST_H
#define _SIM5CONST_H


#define TINY                1e-40                // numerical convenince constant



// sim5 standard unit system is CGS (cm,g,s)
#define grav_radius         1.476716e+05         // gravitational radius GM/c2 of Sun [cm]
#define speed_of_light      2.997925e+10         // speed of light [cm/s]
#define speed_of_light2     8.987554e+20         // square of speed of light [cm^2/s^2]
#define boltzmann_k         1.380650e-16         // Boltzmann constant [erg/K]
#define sb_sigma            5.670400e-05         // Stefan-Boltzmann constant [erg cm-2 s-1 K-4]
#define sigma_thomson       6.652458e-25         // cross-section of Thomson scattering [cm^-2]
#define parsec              3.085680e+18         // parsec [cm]
#define mass_proton         1.672622e-24         // mass of proton [g]
#define mass_electron       9.109382e-28         // mass of electron [g]
#define solar_mass          1.988920e+33         // solar mass [g]
#define grav_const          6.673000e-08         // gravitational constant [cm3 g-1 s-2]
#define planck_h            6.626069e-27         // planck constant [erg.s]

// astrophysical constants
#define Mdot_Edd        2.225475942e+18      // Eddington mass accretion rate [g/s * (M/Msun)]   -- 64*pi*(mass of sun)*G/(speed of light)/(kappa_es); kappa_es=(0.20*(1+X)cm^2/g; X=hydrogen mass fraction (X=1.0)
#define L_Edd           1.257142540e+38      // Eddington luminosity [erg/s * (M/Msun)]  -- 4*pi*(mass of sun)*G*(mass of proton)*(speed of light)/(Thomson sigma)



// constants in SI units (m,kg,s)
#define si_grav_radius      1.476716e+03         // gravitational radius GM/c2 of Sun [m]
#define si_speed_of_light   2.997925e+08         // speed of light [m/s]
#define si_speed_of_light2  8.987554e+16         // square of speed of light [m^2/s^2]
#define si_boltzmann_k      1.380650e-23         // Boltzmann constant [J/K]
#define si_sb_sigma         5.670400e-08         // Stefan-Boltzmann constant [J m-2 s-1 K-4]
#define si_electronvolt     1.602177e-19         // electronvolt [J/elecronvolt]
#define si_parsec           3.085680e+16         // parsec [m]
#define si_grav_const       6.673000e-11         // gravitational constant [m3 kg-1 s-2]
#define si_erg              1.000000e-07         // erg [J]
#define si_solar_mass       1.988920e+30         // solar mass [kg]
#define si_angstrom         1.000000e-10         // angstrom [meters]
#define si_sigma_T          6.652459e-29         // Thomson scattering cross-section for electron [m^2]
#define si_mass_proton      1.672622e-27         // mass of proton [kg]
#define si_mass_electron    9.109382e-31         // mass of electron [kg]
#define si_planck_h         6.626069e-34         // planck constant [J.s]


// constants in geometrical units (G=c=1)
#define gu_sb_sigma         1.562539e-60         // Stefan-Boltzmann constant (sb_sigma/(c^5/G)) [K-4]

// unit conversions
#define erg2kev             6.241507e+08         // erg->keV (erg/(1e3*electronvolt))
#define kev2erg             1.602177e-09         // keV-erg ((1e3*electronvolt)/erg)
#define joule2kev           6.241507e+15         // Joule->keV (1/(1e3*electronvolt))
#define joule2erg           1.000000e+07         // Joule->erg (10^7)
#define erg2joule           1.000000e-07         // erg->Joule (10^-7)
#define kev2joule           1.602177e-16         // keV->Joule (1e3*electronvolt)
#define freq2kev            4.135667e-18         // Hz->keV (planck_h/(1e3*electronvolt))
#define kev2freq            2.417990e+17         // keV->Hz (1e3*electronvolt/planck_h)
#define msq2cmsq            1.000000e+04         // squared centimeters in one squared meter
#define cmsq2msq            1.000000e-04         // squared meters in one squared centimeter
#define kelvin2kev          8.617342e-08         // Kelvin->keV (boltzmann_k/(1e3*electronvolt))
#define kev2kelvin          1.160451e+07         // keV->Kelvin ((1e3*electronvolt)/boltzmann_k)
#define m2cm                1.000000e+02         // meter->centimeter
#define cm2m                1.000000e-02         // centimeter->meter
#define kev2ev              1.000000e+03         // kev->ev 
#define ev2kev              1.000000e-03         // ev->kev 


#endif

