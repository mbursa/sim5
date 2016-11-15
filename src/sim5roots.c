//************************************************************************
//    SIM5 library
//    sim5roots.c - root finding
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute of the Czech Academy of Sciences
//************************************************************************


//! \file sim5roots.c
//! Root finding.
//! 
//! Provides routines for ffinding roots of functions.

// created 2005/07/18

/*
#include <math.h>
#include "sim5utils.h"
*/

//! \cond SKIP
#define MAX_STEPS  500
//! \endcond

long rtbis(double x1, double x2, double xacc, double (*fx)(double), double* result)
//! Root finding.
//! Finds root of a function on an interval. Using bisection method, finds the root of a 
//! function `fx` that is known to lie between `x1` and `x2`. The root, returned as `result`, 
//! will be refined until its accuracy is +/- xacc.
//!
//! @param x1 left boundary of the interval where the root is searched for
//! @param x2 right boundary of the interval where the root is searched for
//! @param xacc accuracy 
//! @param fx function 
//!
//! @result Returns 1 if OK and the root position in `result`, 0 if error.
{
	double dx, f, fmid, xmid, rtb;
	long j=0;

	fmid = (*fx)(x2);
	f    = (*fx)(x1);
	if ((f*fmid) >= 0.0) return(0);//error("rtbis: root is not bracketed");

	if (f < 0.0) {
		rtb = x1;
		dx  = x2-x1;
	}
	else {
		rtb = x2;
		dx  = x1-x2;
	}

	for (j=0; j<MAX_STEPS; j++) {
		dx = dx*0.5;
		xmid = rtb+dx;
		fmid = (*fx)(xmid);
		if (fmid <= 0.0) rtb = xmid;
		if ((fabs(dx) < xacc) || (fmid == 0.0)) break;
	}
	if (j >= MAX_STEPS) {
		error("rtbis: too many steps");
		return(0);
	}

	*result = rtb;
	return(1);
}


#undef MAX_STEPS
