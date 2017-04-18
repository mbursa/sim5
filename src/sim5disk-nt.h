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

#ifndef CUDA

#define DISK_NT_OPTION_LUMINOSITY     1

DEVICEFUNC int disk_nt_setup(double M, double a, double mdot_or_L, double alpha, int options);
DEVICEFUNC void disk_nt_finish();
DEVICEFUNC double disk_nt_r_min();
DEVICEFUNC double disk_nt_flux(double r);
DEVICEFUNC double disk_nt_lum();
DEVICEFUNC double disk_nt_mdot();
DEVICEFUNC double disk_nt_temp(double r);
DEVICEFUNC double disk_nt_sigma(double r);
DEVICEFUNC double disk_nt_ell(double r);
DEVICEFUNC double disk_nt_vr(double r);
DEVICEFUNC double disk_nt_h(double r);
DEVICEFUNC double disk_nt_dhdr(double r);
DEVICEFUNC void disk_nt_dump();

#endif //CUDA

#endif



