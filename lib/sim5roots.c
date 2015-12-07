// $Id: rtbis.c,v 1.3 2005/07/18 12:52:20 bursa Exp $

/*
#include <math.h>
#include "sim5utils.h"
*/
#define MAX_STEPS  500

/*
 * Using bisection, find the root of a function func known to lie
 * between x1 and x2.Theroot, returned as rtbis, will be refined
 * until its accuracy is +/- xacc.
*/

long rtbis(double x1, double x2, double xacc, double (*fx)(double), double* result)
//***************************************************
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
	if (j >= MAX_STEPS) return(0);//error("rtbis: too many steps"); 
	
	*result = rtb;
	return(1);
}


#undef MAX_STEPS
