//************************************************************************
//    sim5disk-nt.h - model for N-T disk radial structure
//------------------------------------------------------------------------
//    Date   : 12.10.2014
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2014 Michal Bursa
//************************************************************************


#ifndef _SIM5DISKNT_H
#define _SIM5DISKNT_H

#define DISK_NT_OPTION_LUMINOSITY      1


int disk_nt_setup(double M, double a, double mdot_L, double alpha, int _options);
void disk_nt_finish();
double disk_nt_r_min();
double disk_nt_flux(double r);
double disk_nt_lum();
double disk_nt_mdot();
double disk_nt_temp(double r);
double disk_nt_sigma(double r);
double disk_nt_ell(double r);
double disk_nt_vr(double r);
double disk_nt_h(double r);
double disk_nt_dhdr(double r);
void disk_nt_dump();


#endif



