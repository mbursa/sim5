//************************************************************************
//    sim5rootfinding.h - data file
//------------------------------------------------------------------------
//    Date   : 12.10.2004
//    Author : Michal Bursa
//    e-mail : bursa@sirrah.troja.mff.cuni.cz
//------------------------------------------------------------------------
//    (C) 2004 Michal Bursa
//************************************************************************
// $Id: sim4datafile.h,v 1.7 2005/12/08 10:36:12 bursa Exp $

#ifndef _SIM5ROOTFINDING_H
#define _SIM5ROOTFINDING_H

#ifdef __cplusplus
extern "C" {
#endif

long rtbis(double x1, double x2, double xacc, double (*fx)(double), double* result);

#ifdef __cplusplus
}
#endif

#endif
