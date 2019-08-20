//************************************************************************
//    sim5elliptic.h - data file
//------------------------------------------------------------------------
//    Date   : 12.10.2004
//    Author : Michal Bursa
//    e-mail : bursa@sirrah.troja.mff.cuni.cz
//------------------------------------------------------------------------
//    (C) 2004 Michal Bursa
//************************************************************************
// $Id: sim4datafile.h,v 1.7 2005/12/08 10:36:12 bursa Exp $

#ifndef _SIM5ELLIPTIC_H
#define _SIM5ELLIPTIC_H


#ifdef __cplusplus
extern "C" {
#endif

DEVICEFUNC double rf(double x, double y, double z);
DEVICEFUNC double rd(double x, double y, double z);
DEVICEFUNC double rc(double x, double y);
DEVICEFUNC double rj(double x, double y, double z, double p);
DEVICEFUNC INLINE double elliptic_k(double m);
DEVICEFUNC double elliptic_e_sin(double sin_phi, double m);
DEVICEFUNC double elliptic_e_cos(double cos_phi, double m);
DEVICEFUNC double elliptic_f(double phi, double m);
DEVICEFUNC double elliptic_f_cos(double cos_phi, double m);
DEVICEFUNC double elliptic_f_sin(double sin_phi, double m);
DEVICEFUNC double elliptic_pi_complete(double n, double m);
DEVICEFUNC double elliptic_pi_sin(double sin_phi, double n, double m);
DEVICEFUNC sim5complex elliptic_pi(double phi, double n, double m);
DEVICEFUNC INLINE double jacobi_isn(double z, double m);
DEVICEFUNC INLINE double jacobi_icn(double z, double m);
DEVICEFUNC INLINE double jacobi_itn(double z, double m);
DEVICEFUNC void  jacobi_sncndn(double u, double m, double *sn, double *cn, double *dn);
DEVICEFUNC INLINE double jacobi_sn(double u, double m);
DEVICEFUNC INLINE double jacobi_cn(double u, double m);
DEVICEFUNC INLINE double jacobi_dn(double u, double m);

DEVICEFUNC double integral_R_r0_re(double a, double b, double c, double d, double X);
DEVICEFUNC double integral_R_r0_cc(double a, double b, sim5complex c, double X);
DEVICEFUNC double integral_R_r0_re_inf(double a, double b, double c, double d);
DEVICEFUNC double integral_R_r0_cc_inf(double a, double b, sim5complex c);
DEVICEFUNC double integral_R_r1_re(double a, double b, double c, double d, double X);
DEVICEFUNC double integral_R_r1_cc(double a, double b, sim5complex c, double X1, double X2);
DEVICEFUNC double integral_R_r2_re(double a, double b, double c, double d, double X);
DEVICEFUNC double integral_R_r2_cc(double a, double b, sim5complex c, double X1, double X2);
DEVICEFUNC double integral_R_rp_re(double a, double b, double c, double d, double p, double X);
DEVICEFUNC double integral_R_rp_cc2(double a, double b, sim5complex c, double p, double X1, double X2);
DEVICEFUNC double integral_R_rp_re_inf(double a, double b, double c, double d, double p);
DEVICEFUNC double integral_R_rp_cc2_inf(double a, double b, sim5complex c, double p, double X1);
DEVICEFUNC double integral_T_m0(double a2, double b2, double X);
DEVICEFUNC double integral_T_m2(double a2, double b2, double X);
DEVICEFUNC double integral_T_mp(double a2, double b2, double p, double X);

#ifdef __cplusplus
}
#endif

#endif
