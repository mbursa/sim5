/*
#include "sim5config.h"
#ifndef CUDA
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "sim5math.h"
#include "sim5utils.h"
#include "sim5elliptic.h"
#endif
*/

// rf(x,y,z) computes Calson's elliptic integral of the first kind, R_F(x,y,z).
// x, y, and z must be nonnegative, and at most one can be zero. See "Numerical
// Recipes in C" by W. H. Press et al. (Chapter 6)
DEVICEFUNC
static double rf(double x, double y, double z)
{
	const double ERRTOL=0.0003, rfTINY=1.5e-308, rfBIG=3.0e307, THIRD=1.0/3.0;
	const double C1=1.0/24.0, C2=0.1, C3=3.0/44.0, C4=1.0/14.0;
	double alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt;

	if ((min(min(x,y),z) < 0.0) || (min(min(x+y,x+z),y+z) < rfTINY) || (max(max(x,y),z) > rfBIG)) {
        #ifndef CUDA
		warning("%e/%e/%e\n", x,y,z);
		error("invalid arguments in rf");
        #endif
	}

	xt=x;
	yt=y;
	zt=z;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=THIRD*(xt+yt+zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (max(max(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	e2=delx*dely-delz*delz;
	e3=delx*dely*delz;
	return (1.0+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave);
}


// rd(x,y,z) computes Calson's elliptic integral of the second kind, R_D(x,y,z).
// x, y, and z must be nonnegative, and at most one can be zero. z must be positive. 
// See "Numerical Recipes in C" by W. H. Press et al. (Chapter 6)
DEVICEFUNC
static double rd(double x, double y, double z) {
	const double ERRTOL=0.0003, rdTINY=1.5e-308, rdBIG=3.0e307;
	const double C1=3.0/14.0, C2=1.0/6.0, C3=9.0/22.0, C4=3.0/26.0, C5=0.25*C3, C6=1.5*C4;
	double alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt;

	if ((min(x,y) < 0.0) || (min(x+y,z) < rdTINY) || (max(max(x,y),z) > rdBIG)) {
        #ifndef CUDA
		fprintf(stderr, "%e/%e/%e\n", x,y,z);
		error("invalid arguments in rd");
        #endif
	}

	xt=x;
	yt=y;
	zt=z;
	sum=0.0;
	fac=1.0;
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		sum += fac/(sqrtz*(zt+alamb));
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		ave=0.2*(xt+yt+3.0*zt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
	} while (max(max(fabs(delx),fabs(dely)),fabs(delz)) > ERRTOL);
	ea=delx*dely;
	eb=delz*delz;
	ec=ea-eb;
	ed=ea-6.0*eb;
	ee=ed+ec+ec;
	return 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*delz*ee)
		+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave));
}


// Computes Carlson’s degenerate elliptic integral, RC(x,y) 
// x must be nonnegative and y must be nonzero. 
// If y < 0, the Cauchy principal value is returned.
DEVICEFUNC
static double rc(double x, double y) {
	const double ERRTOL=0.0003, rcTINY=1.5e-308, rcBIG=3.0e37;
	const double THIRD=1.0/3.0, C1=0.3, C2=1.0/7.0, C3=0.375, C4=9.0/22.0;
	const double COMP1=2.236/sqrt(rcTINY),COMP2=sqr(rcTINY*rcBIG)/25.0;
	
	double alamb,ave,s,w,xt,yt;

	if ((x < 0.0) || (y == 0.0) || ((x+fabs(y)) < rcTINY) || ((x+fabs(y)) > rcBIG) ||
	  ((y<-COMP1) && (x > 0.0) && (x < COMP2))) {
        #ifndef CUDA
		fprintf(stderr, "%e/%e\n", x,y);
		error("invalid arguments in rc");
        #endif
	}
	
	if (y > 0.0) {
		xt=x;
		yt=y;
		w=1.0;
	} else {
		xt=x-y;
		yt= -y;
		w=sqrt(x)/sqrt(xt);
	}
	do {
		alamb=2.0*sqrt(xt)*sqrt(yt)+yt;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		ave=THIRD*(xt+yt+yt);
		s=(yt-ave)/ave;
	} while (fabs(s) > ERRTOL);
	return w*(1.0+s*s*(C1+s*(C2+s*(C3+s*C4))))/sqrt(ave);
}


// rj(x,y,z) computes Calson's elliptic integral of the third kind, R_J(x,y,z).
// x, y, and z must be nonnegative, and at most one can be zero. p must be nonzero.
// If p < 0, the Cauchy principal value is returned. 
// See "Numerical Recipes in C" by W. H. Press et al. (Chapter 6)
DEVICEFUNC
static double rj(double x, double y, double z, double p) {
	const double ERRTOL=0.0003, rjTINY=pow(5.0*1.5e-308,1./3.), rjBIG=0.3*pow(0.2*3.0e307,1./3.);
	const double C1=3.0/14.0, C2=1.0/3.0, C3=3.0/22.0, C4=3.0/26.0, 
		C5=0.75*C3, C6=1.5*C4, C7=0.5*C2, C8=C3+C3;
	double a,alamb,alpha,ans,ave,b,beta,delp,delx,dely,delz,ea,eb,ec,ed,ee,
		fac,pt,rcx,rho,sqrtx,sqrty,sqrtz,sum,tau,xt,yt,zt;
	if ((min(min(x,y),z) < 0.0) || (min(min(x+y,x+z),min(y+z,fabs(p))) < rjTINY) ||
	  (max(max(x,y),max(z,fabs(p))) > rjBIG)) {
        #ifndef CUDA
		fprintf(stderr, "WRN: invalid arguments in rj (%e/%e/%e/%e\n)", x,y,z,p);
		//error("invalid arguments in rj");
        #endif
        return 0.0;
	}

    a=b=rcx=0.0;
	sum=0.0;
	fac=1.0;
	if (p > 0.0) {
		xt=x;
		yt=y;
		zt=z;
		pt=p;
	} else {
		xt=min(min(x,y),z);
		zt=max(max(x,y),z);
		yt=x+y+z-xt-zt;
		a=1.0/(yt-p);
		b=a*(zt-yt)*(yt-xt);
		pt=yt+b;
		rho=xt*zt/yt;
		tau=p*pt/yt;
		rcx=rc(rho,tau);
	}
	do {
		sqrtx=sqrt(xt);
		sqrty=sqrt(yt);
		sqrtz=sqrt(zt);
		alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
		alpha=sqr(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
		beta=pt*sqr(pt+alamb);
		sum += fac*rc(alpha,beta);
		fac=0.25*fac;
		xt=0.25*(xt+alamb);
		yt=0.25*(yt+alamb);
		zt=0.25*(zt+alamb);
		pt=0.25*(pt+alamb);
		ave=0.2*(xt+yt+zt+pt+pt);
		delx=(ave-xt)/ave;
		dely=(ave-yt)/ave;
		delz=(ave-zt)/ave;
		delp=(ave-pt)/ave;
	} while (max(max(fabs(delx),fabs(dely)),max(fabs(delz),fabs(delp))) > ERRTOL);
	ea=delx*(dely+delz)+dely*delz;
	eb=delx*dely*delz;
	ec=delp*delp;
	ed=ea-3.0*ec;
	ee=eb+2.0*delp*(ea-ec);
	ans=3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))
		+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
	if (p <= 0.0) ans=a*(b*ans+3.0*(rcx-rf(xt,yt,zt)));
	return ans;
}

//------------------------------------------------------------------------------------


// elliptic_k(m) = K(m) = F(pi/2,m), the complete elliptic integral of 
// the first kind. Note, here m=k^2, k is used in P. F. Byrd & M. D. 
// Friedman, Handbook of Elliptic Integrals for Engineers and Physicists. [This 
// subroutine has been tested with Mathematica 4.2. Note, m is also used in 
// Mathematica and the paper of Cadez et al.] The subroutine is applicable only 
// if 0 <= m <1.
DEVICEFUNC INLINE
double elliptic_k(double m)
{
  if (m==1.0) m = 0.99999999;
  #ifndef CUDA
  if (m>1.0) error("Error in routine elliptic_k(m): m >= 1.0 (%e)", m);
  #endif
  return rf(0,1.0-m,1.0);
}


//-----------------------------------------------------------------------------------

// elliptic_f(phi,m) = F(phi,m), the complete elliptic integral of 
// the first kind. Note, here m=k^2, k is used in P. F. Byrd & M. D. 
// Friedman, Handbook of Elliptic Integrals for Engineers and Physicists. [This 
// subroutine has been tested with Mathematica 4.2. Note, m is also used in 
// Mathematica and the paper of Cadez et al.] The subroutine is applicable only 
// if 0 <= m <1.
DEVICEFUNC
double elliptic_f(double phi, double m)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_f(m): m >= 1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (phi==0.0) return 0.0;

    int k = 0;
    while (fabs(phi) > M_PI/2.) { (phi>0)?k++:k--; phi += (phi>0)?-M_PI:+M_PI;  }

	double s2 = pow(sin(phi), 2);
	double ell_f = (phi>0?+1:-1) * sqrt(s2)*rf(1-s2, 1.0-s2*m, 1.0);
    if (k!=0) ell_f += 2.*k*elliptic_k(m);
    return ell_f;
}

DEVICEFUNC
double elliptic_f_cos(double cos_phi, double m)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_f_cos(m): m >= 1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (cos_phi==1.0) return 0.0;

	double X = 0.0;
	if (cos_phi<0.0) {
		cos_phi = - cos_phi;
		X = 2.0 * rf(0.0, 1.0-m, 1.0);
	} 
		
	double s2 = 1.0-sqr(cos_phi);
	return X + ((X==0.0)?(+1):(-1))*sqrt(s2)*rf(1.0-s2, 1.0-s2*m, 1.0);
}

DEVICEFUNC
double elliptic_f_sin(double sin_phi, double m)
{
    #ifndef CUDA
	if ((sin_phi<0.0)||(sin_phi>1.0)) error("elliptic_f_sin: invalid input (sin_phi<0)");
	if (m>1.0) error("Error in routine elliptic_f_sin(m): m >= 1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (sin_phi==0.0) return 0.0;
	double s2 = sqr(sin_phi);
	return sin_phi*rf(1.-s2, 1.0-s2*m, 1.0);
}


//-----------------------------------------------------------------------------------

// Legendre elliptic integral of the second kind E(phi,m), evaluated using Carlson’s
// functions RD and RF . 
// The argument ranges are 0 <= phi <= pi/2; 0 <= sqrt(m)*sin(phi) <= 1


DEVICEFUNC
double elliptic_e(double phi, double m)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_e(m): m >= 1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (phi==0.0) return 0.0;

    #ifndef CUDA
	if ((phi<0.0)||(phi>M_PI)) error("elliptic_e: invalid input");
    #endif
	double X = 0.0;
	if (phi>M_PI/2.) {
		phi = M_PI - phi;
		X = 2.0 * ( rf(0.0,1.0-m,1.0) - m*rd(0.0,1.0-m,1.0)/3.0 );
	} 
		
	double s  = sin(phi);
	double c2 = sqr(cos(phi));
	double q  = 1.0 - s*s*m;
	return X + ((X==0.0)?(+1):(-1))*s*( rf(c2,q,1.0) - sqr(s*sqrt(m))*rd(c2,q,1.0)/3.0 );
}

DEVICEFUNC
double elliptic_e_cos(double cos_phi, double m)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_e_cos: m >= 1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (cos_phi==1.0) return 0.0;

	double X = 0.0;
	if (cos_phi<0.0) {
		cos_phi = -cos_phi;
		X = 2.0 * ( rf(0.0,1.0-m,1.0) - m*rd(0.0,1.0-m,1.0)/3.0 );
	} 
		
	double c2 = sqr(cos_phi);
	double s  = sqrt(1.0-c2);
	double q  = 1.0 - m + c2*m;
	return X + ((X==0.0)?(+1):(-1))*s*( rf(c2,q,1.0) - sqr(s*sqrt(m))*rd(c2,q,1.0)/3.0);
}

DEVICEFUNC
double elliptic_e_sin(double sin_phi, double m)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_e_sin: m >= 1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (sin_phi==0.0) return 0.0;

    #ifndef CUDA
	if ((sin_phi<0.0)||(sin_phi>1.0)) error("elliptic_e_sin: invalid input (sin_phi<0)");
    #endif
	double s2 = sin_phi*sin_phi;
	double c2 = 1.0 - s2;
	double q  = 1.0 - s2*m;
	return sin_phi*( rf(c2,q,1.0) - sqr(sin_phi*sqrt(m))*rd(c2,q,1.0)/3.0);
}

//-----------------------------------------------------------------------------------


// Legendre elliptic integral of the third kind PI(phi,n,m), evaluated using Carlson’s
// functions RJ and RF . 
// note: the sign convetion of n is the same as in Mathematica
// The argument ranges are 0 <= phi <= pi/2; 0 <= sqrt(m)*sin(phi) <= 1
// complete integral: phi=pi/2
DEVICEFUNC
double elliptic_pi_complete(double n, double m)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_pi_complete: m>1.0");
	if (n>1.0) error("Error in routine elliptic_pi_complete: n>1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (n==1.0) n = 0.99999999;

	double q  = 1.0-m;
	return rf(0.0,q,1.0) + n*rj(0.0,q,1.0,1.0-n)/3.0;
}


DEVICEFUNC
//double complex elliptic_pi(double phi, double n, double m)
sim5complex elliptic_pi(double phi, double n, double m)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_pi: m>1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (phi==0.0) return makeComplex(0.0,0.0);
//	if ((n*pow(sin(phi),2.))>=1.0) error("Error in routine elliptic_pi: n>=1/sin(phi)^2 (%e/%e/%e)", phi,n,m);
	
	int hyperbolicIII = (n > 1.0); 
    double p=0,dn=0,la=0;
    if (hyperbolicIII) {
        p  = sqrt((n-1.)*(1.-m/n));
        dn = sqrt(1. - m*sin(phi)*sin(phi));
        la = (dn+p*tan(phi)) / (dn-p*tan(phi));
        n = m/n;
    }

    int k = 0;
    while (fabs(phi) > M_PI/2.) { (phi>0)?k++:k--; phi += (phi>0)?-M_PI:+M_PI;  }

    double s = sin(phi);
	double c2  = 1.0 - s*s;
	double q   = 1.0 - s*s*m;
	double ns2 = -n*s*s;

    #ifndef CUDA
    if (rj(c2,q,1.0,1.0+ns2)==0.0) error("Error in routine elliptic_pi: rj==0 (%e/%e/%e)", phi,n,m);
    #endif

    //double complex ell_pi;    
    sim5complex ell_pi = makeComplex(s*(rf(c2,q,1.0) - ns2*rj(c2,q,1.0,1.0+ns2)/3.0), 0.0);
    if ((k!=0) && (!hyperbolicIII)) ell_pi += 2.*k*elliptic_pi_complete(n,m);

    //if (hyperbolicIII) ell_pi = -ell_pi + elliptic_f(phi,m) + log(fabs(la))/(2*p) + ((la<0)?I*M_PI/(2*p):0.0);
    if (hyperbolicIII) ell_pi = -ell_pi + elliptic_f(phi,m) + log(fabs(la))/(2*p) + ((la<0)?ComplexI*M_PI/(2*p):nullComplex());

    return ell_pi;
}


DEVICEFUNC
double elliptic_pi_cos(double cos_phi, double n, double m)
{
	//double s2 = 1.0-cos_phi*cos_phi;

    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_pi_cos: m>1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (cos_phi==1.0) return 0.0;
    //if (fabs(n*s2)>=1.0) error("Error in routine elliptic_pi_cos: n>=1/sin(phi)^2 (%e/%e/%e)", cos_phi,n,m);

	double X = 0.0;
	if (cos_phi<0.0) {
		cos_phi = -cos_phi;
		X = 2.0 * ((rf(0.0,1.0-m,1.0) + n*rj(0.0,1.0-m,1.0,1.0-n)/3.0));
	} 
		
	double c2  = sqr(cos_phi);
	double s   = sqrt(1.0-c2);
	double ns2 = -n*(1.0-c2);
	double q  = 1.0 - (1.0-c2)*m;
	return X + ((X==0.0)?(+1):(-1))*s*(rf(c2,q,1.0) - ns2*rj(c2,q,1.0,1.0+ns2)/3.0);
}


DEVICEFUNC
double elliptic_pi_sin(double sin_phi, double n, double m)
{
	double s2 = sin_phi*sin_phi;

    #ifndef CUDA
	if (m>1.0) error("Error in routine elliptic_pi_sin: m>1.0");
    #endif
	if (m==1.0) m = 0.99999999;
	if (sin_phi==0.0) return 0.0;
	//if (fabs(n*s2)>=1.0) error("Error in routine elliptic_pi_sin: n>=1/sin(phi)^2 (%e/%e/%e)", sin_phi,n,m);

    #ifndef CUDA
	if ((sin_phi<0.0)||(sin_phi>1.0)) error("elliptic_pi_sin: invalid input (sin_phi<0)");
    #endif
	double c2  = 1.0-s2;
	double ns2 = -n*s2;
	double q  = 1.0 - s2*m;
	return sin_phi*(rf(c2,q,1.0) - ns2*rj(c2,q,1.0,1.0+ns2)/3.0);
}


// jacobi_isn(y,emmc) = sn^(-1) (y, emmc), the inverse Jacobian elliptic function 
// sn^(-1). Note emmc=m = k^2. [This subroutine has been tested with Mathematica.]
// Applicable only if 0 <= emmc <1 and -1 < y <1.
DEVICEFUNC INLINE
double jacobi_isn(double y, double emmc)
{
    return y*rf(1-y*y,1.0-emmc*y*y,1.0);
}


DEVICEFUNC INLINE
double jacobi_icn1(double z, double m)
{
    if (fabs(z-1.0)<1e-8) return 0.0;
    if (fabs(m-0.0)<1e-8) return acos(z);
    #ifndef CUDA
    if (fabs(z)>1.0) error("jacobi_icn1: z>1 (%.10e)", z);
    if (m>=1.0) error("jacobi_icn: x<0 and m>=1.0");
    #endif
    double z2 = z*z;
    return sqrt(1-z2) * rf(z2, 1.0 - m*(1.-z2), 1.0);
}


// jacobi_icn(y,emmc) = cn^(-1) (y, emmc), the inverse Jacobian elliptic function 
// cn^(-1). Note emmc=m = k^2. [This subroutine has been tested with Mathematica.]
// Applicable only if 0 <= m <1 and -1 <= z <= 1.
DEVICEFUNC INLINE
double jacobi_icn(double z, double m)
{
    if (z==0.0) {
        return elliptic_k(m);
    } else
    if (z>0.0) {
        return jacobi_icn1(z,m);
    } else {
        return 2.*elliptic_k(m)-jacobi_icn1(-z,m);
    }
}


// sncndn(uu,emmc,sn,cn,dn) returns the Jacobian elliptic functions sn(uu,1-emmc), 
// cn(uu,1-emmc), and dn(uu,1-emmc). See "Numerical Recipes in C" by W. H. Press 
// et al. (Chapter 6).
DEVICEFUNC
void sncndn(double uu, double emmc, double *sn, double *cn, double *dn)
{
    #ifndef CUDA
	if (emmc>1.0) error("Error in routine sncndn: m>1.0");
	if (uu>2.*elliptic_k(emmc)) error("sncndn: invalid input (u>2K(m))");
    #endif
	if (emmc==1.0) emmc = 0.99999999;

	const double CA=1.0e-8;
	int bo;
	int i,ii,l;
	double a,b,c,d,emc,u;
	double em[13],en[13];

    d=1.0;
	emc=1.0-emmc;
	u=uu;
	if (emc != 0.0) {
		bo=(emc < 0.0);
		if (bo) {
			d=1.0-emc;
			emc /= -1.0/d;
			u *= (d=sqrt(d));
		}
		a=1.0;
		*dn=1.0;
		for (i=0;i<13;i++) {
			l=i;
			em[i]=a;
			en[i]=(emc=sqrt(emc));
			c=0.5*(a+emc);
			if (fabs(a-emc) <= CA*a) break;
			emc *= a;
			a=c;
		}
		u *= c;
		*sn=sin(u);
		*cn=cos(u);
		if (*sn != 0.0) {
			a=(*cn)/(*sn);
			c *= a;
			for (ii=l;ii>=0;ii--) {
				b=em[ii];
				a *= c;
				c *= *dn;
				*dn=(en[ii]+a)/(b+a);
				a=c/b;
			}
			a=1.0/sqrt(c*c+1.0);
			*sn=((*sn) >= 0.0 ? a : -a);
			*cn=c*(*sn);
		}
		if (bo) {
			a=*dn;
			*dn=*cn;
			*cn=a;
			*sn /= d;
		}
	} else {
		*cn=1.0/cosh(u);
		*dn=*cn;
		*sn=tanh(u);
	}
}

// jacobi_sn(uu,emmc) = sn(uu,emmc), the Jacobian elliptic function sn. Note, here 
// emmc= m = k^2, k is used in P. F. Byrd & M. D. Friedman, Handbook of Elliptic 
// Integrals for Engineers and Physicists. The subroutine has been tested with 
// Mathematica 4.2. 
DEVICEFUNC INLINE
double jacobi_sn(double uu, double emmc)
{
  double snx, cnx, dnx;

  sncndn(uu,emmc,&snx,&cnx,&dnx);
  return snx;
}

// jacobi_cn(uu,emmc) = cn(uu,emmc), the Jacobian elliptic function cn. Note, here 
// emmc= m = k^2. The subroutine has been tested with Mathematica 4.2.
DEVICEFUNC INLINE
double jacobi_cn(double uu, double emmc)
{
  double snx, cnx, dnx;

  sncndn(uu,emmc,&snx,&cnx,&dnx);
  return cnx;
}


// Jacobian elliptic function dn
DEVICEFUNC INLINE
double jacobi_dn(double u, double m)
{
  double sn, cn, dn;
  sncndn(u,m,&sn,&cn,&dn);
  return dn;
}



//---------------------------------------------------------------------------

DEVICEFUNC INLINE
double integral_C0(double u, double m)
// int du
// Eq. 312.00 (Byrd & Friedman)
{
	return u;
}


DEVICEFUNC INLINE
double integral_C1(double u, double m)
// int cn(u) du
// Eq. 312.01 (Byrd & Friedman)
{
	double sn, cn, dn;
	sncndn(u,m,&sn,&cn,&dn);
	return acos(dn)/sqrt(m);
}


DEVICEFUNC INLINE
double integral_C2(double u, double m)
// int cn(u)^2 du
// Eq. 312.02 (Byrd & Friedman)
{
	double sn, cn, dn;
	sncndn(u,m,&sn,&cn,&dn);
	return 1./m*(elliptic_e_cos(cn,m) - (1.-m)*u);
}


DEVICEFUNC INLINE
double integral_C2_cos(double cn_u, double m)
// int cn(u)^2 du
// Eq. 312.02 (Byrd & Friedman)
{
	return 1./m*(elliptic_e_cos(cn_u,m) - (1.-m)*elliptic_f_cos(cn_u,m));
}


DEVICEFUNC
double integral_Z1(double a, double b, double u, double m)
// int (1-b*sn(u,m)^2)/(1-a*sn(u,m)^2) du
// a,b > 0
// Eq. 340.01 (Byrd & Friedman)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine integral_Z1 (m>1.0)");
//	if ((a<0)||(b<0)) error("integral_Z1: invalid input (a|b<0)");
	if (u>2.*elliptic_k(m)) error("integral_Z1: invalid input (u>2K(m))");
    #endif
	double sn, cn, dn;
	sncndn(u, m, &sn, &cn, &dn);
	return 1./a*((a-b)*elliptic_pi_cos(cn,a,m) + b*u);
}


DEVICEFUNC
double integral_Z2(double a, double b, double u, double m)
// int (1-b*sn(u,m)^2)^2/(1-a*sn(u,m)^2)^2 du
// a,b > 0
// Eq. 340.02 (Byrd & Friedman)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine integral_Z2: m>1.0");
//	if ((a<0)||(b<0)) error("integral_Z2: invalid input (a|b<0)");
	if (u>2.*elliptic_k(m)) error("integral_Z2: invalid input (u>2K(m))");
    #endif
	double sn, cn, dn;
	sncndn(u, m, &sn, &cn, &dn);
	//fprintf(stderr,"cn=%e  a=%e  m=%e\n", cn,a,m);
	double V1 = elliptic_pi_cos(cn,a,m);
	double V2 = 0.5/((a-1.)*(m-a)) * ( 
	              a*elliptic_e_cos(cn,m) + (m-a)*u + 
	              (2.*a*m+2.*a-a*a-3.*m)*V1 - 
	              (a*a*sn*cn*dn)/(1.-a*sn*sn)
               );
        double ab = a-b;
	return 1./sqr(a)*( sqr(b)*u + 2.*b*ab*V1 + ab*ab*V2 );
}


DEVICEFUNC INLINE
double integral_Rm1(double a, double u, double m)
// int (1+a*cn(u)) du
// Eq. 341.00 (Byrd & Friedman)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine integral_Rm1: m>1.0");
	//if (u>2.*elliptic_k(m)) error("integral_Rm1: invalid input (u>2K(m))");
    #endif

	return u + a/sqrt(m)*acos(jacobi_dn(u,m));
}


DEVICEFUNC
double integral_Rm2(double a, double u, double m)
// int (1+a*cn(u))^2 du
// Eq. 341.01 (Byrd & Friedman)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine integral_Rm2: m>1.0");
	//if (u>2.*elliptic_k(m)) error("integral_Rm2: invalid input (u>2K(m))");
    #endif
    double a2 = sqr(a);
	double sn, cn, dn;
	sncndn(u, m, &sn, &cn, &dn);
	return 1/m*( (m-a2*(1.-m))*u + a2*elliptic_e_cos(cn,m) + 2*a*sqrt(m)*acos(dn) );
}


DEVICEFUNC INLINE
double integral_R0(double u, double m)
{
	return u;
}


DEVICEFUNC
double integral_R1(double a, double u, double m)
// int 1/(1+a*cn(u)) du
// a != 1.0
// Eq. 341.03, 361.54 (Byrd & Friedman)
// Note a problem with 361.54, where the f1 expression for a2(a2-1)>m case is wrong
// and also the sign of a*f1 seem to be oposite (we have: ellpi + a*f1) 
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine integral_R1: m>1.0");
	//if (u>2.*elliptic_k(m)) error("integral_R1: invalid input (u>2K(m))");
    #endif
	
	double a2 = sqr(a);
	double n = a2/(a2-1.);
	double sn, cn, dn;
	sncndn(u, m, &sn, &cn, &dn);
	double mma = (m + (1.-m)*a2) / (1.-a2);
	sim5complex f1    = (fabs(mma)>1e-5) ? csqrt(makeComplex(1./mma,0.0))*catan(csqrt(makeComplex(mma,0.0))*sn/dn) : makeComplex(sn/dn,0.0);
	sim5complex ellpi = elliptic_pi_cos(cn, n, m);
	sim5complex res   = 1./(1.-a2) * (ellpi + a*f1);
	/* 
    fprintf(stderr, "R1 - a   = %f\n",a);
	fprintf(stderr, "R1 - u   = %f\n",u);
	fprintf(stderr, "R1 - m   = %f\n",m);
    fprintf(stderr, "R1 - f1  = %e + %e I\n", creal(f1), cimag(f1));
    fprintf(stderr, "R1 - af1  = %e + %e I\n", a*creal(f1), a*cimag(f1));
    fprintf(stderr, "R1 - epi = %e + %e I\n", creal(ellpi), cimag(ellpi));
    fprintf(stderr, "R1 - res = %e + %e I\n", creal(res), cimag(res));
    fprintf(stderr, "R1-tot = %e + %e I\n", -1.114784e-01-creal(res), cimag(res));
     */
	return creal(res);
}


DEVICEFUNC
double integral_R2(double a, double u, double m)
// int 1/(1+a*cn(u))^2 du
// a != 1.0
// Eq. 341.04 (Byrd & Friedman)
{
    #ifndef CUDA
	if (m>1.0) error("Error in routine integral_R2: m>1.0");
	//if (u>2.*elliptic_k(m)) error("integral_R2: invalid input (u>2K(m))");
    #endif
	
    double a2  = sqr(a);
	double mma = (m + (1.-m)*a2);
	double sn, cn, dn;
	sncndn(u, m, &sn, &cn, &dn);
	
	return 1/(a2-1.)/mma * ( 
	    (a2*(2.*m-1.)-2.*m)*integral_R1(a,u,m) + 
	    2.*m*integral_Rm1(a,u,m) - 
	    m*integral_Rm2(a,u,m) + 
	    a*a2*sn*dn/(1.+a*cn) 
    );
}







DEVICEFUNC
double integral_R_r0_re(double a, double b, double c, double d, double X)
// int_a^X 1/sqrt((x-a)(x-b)(x-c)(x-d))
// X > a > b > c > d
// Eq. 258.00 (Byrd & Friedman)
{
    #ifndef CUDA
	if ((a<=b)||(b<=c)||(c<=d)) error("integral_R_r0_re: invalid input");
    #endif
	double m4 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt(((b-d)*(X-a))/((a-d)*(X-b))); 
	//return 2.0/sqrt((a-c)*(b-d)) * elliptic_f_sin(sn, m);  // ... should be same as bellow
	return 2.0/sqrt((a-c)*(b-d)) * jacobi_isn(sn, m4);       // but this is faster to compute
}


DEVICEFUNC
double integral_R_r0_re_inf(double a, double b, double c, double d)
// int_a^infinity 1/sqrt((x-a)(x-b)(x-c)(x-d))
// X > a > b > c > d
// Eq. 258.00 (Byrd & Friedman)
{
    #ifndef CUDA
	if ((a<=b)||(b<=c)||(c<=d)) error("integral_R_r0_re: invalid input");
    #endif
	double m4 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt((b-d)/(a-d)); 
	//return 2.0/sqrt((a-c)*(b-d)) * elliptic_f_sin(sn, m);
	return 2.0/sqrt((a-c)*(b-d)) * jacobi_isn(sn, m4);
}


DEVICEFUNC
double integral_R_r0_cc(double a, double b, sim5complex c, double X)
// int_a^X 1/sqrt((x-a)(x-b)(x-c)(x-d))
// X > a > b;  c = u + i*v;  d = c* 
// Eq. 260.00 (Byrd & Friedman)
{
	double u  = creal(c); 
	double v2 = sqr(cimag(c)); 
	double A  = sqrt(sqr(a-u) + v2); 
	double B  = sqrt(sqr(b-u) + v2); 
	double m2 = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double cn = (X*(A-B)+a*B-b*A) / (X*(A+B)-a*B-b*A);
	//return 1./sqrt(A*B) * elliptic_f_cos(cn, m);
	return 1./sqrt(A*B) * jacobi_icn(cn, m2);

}


DEVICEFUNC
double integral_R_r0_cc_inf(double a, double b, sim5complex c)
// int_a^infinity 1/sqrt((x-a)(x-b)(x-c)(x-d))
// X > a > b;  c = u + i*v;  d = c* 
// Eq. 260.00 (Byrd & Friedman)
{
	double u  = creal(c); 
	double v2 = sqr(cimag(c)); 
	double A  = sqrt(sqr(a-u) + v2); 
	double B  = sqrt(sqr(b-u) + v2); 
	double m2 = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double cn = (A-B)/(A+B);
	//return 1./sqrt(A*B) * elliptic_f_cos(cn, m);
	return 1./sqrt(A*B) * jacobi_icn(cn, m2);
}


DEVICEFUNC
double integral_R_r1_re(double a, double b, double c, double d, double X)
// int_a^X x/sqrt((x-a)(x-b)(x-c)(x-d))
// X > a > b > c > d
// Eq. 258.11 (Byrd & Friedman)
{
	double m2 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt(((b-d)*(X-a))/((a-d)*(X-b))); 
	double u  = jacobi_isn(sn, m2);//elliptic_f_sin(sn,m);
	double a2 = (a-d)/(b-d);
	double b2 = ((a-d)*b)/(a*(b-d)); 
	double Z  = integral_Z1(a2,b2,u,m2) - integral_Z1(a2,b2,0,m2);
	return a*2.0/sqrt((a-c)*(b-d)) * Z;
}


DEVICEFUNC
//double integral_R_r1_cc(double a, double b, double complex c, double X1, double X2)
double integral_R_r1_cc(double a, double b, sim5complex c, double X1, double X2)
// int_a^X x/sqrt((x-a)(x-b)(x-c)(x-c*))
// X > a > b;  c = u + i*v;
// Eq. 260.03 (Byrd & Friedman)
{
	double u  = creal(c); 
	double v2 = sqr(cimag(c)); 
	double A  = sqrt(sqr(a-u) + v2); 
	double B  = sqrt(sqr(b-u) + v2); 
	double m  = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double g  = 1./sqrt(A*B);
	double alpha1 = (B*a+b*A)/(B*a-b*A);
	double alpha2 = (B+A)/(B-A);

	double u1 = elliptic_f_cos((X1*(A-B)+a*B-b*A)/(X1*(A+B)-a*B-b*A),m);
	double u2 = elliptic_f_cos((X2*(A-B)+a*B-b*A)/(X2*(A+B)-a*B-b*A),m);

	double t0 = alpha1 * (integral_R0(u2,m) - integral_R0(u1,m));
	double t1 = (alpha2-alpha1) * (integral_R1(alpha2,u2,m) - integral_R1(alpha2,u1,m));

	return (B*a-b*A)/(B+A)*g * (t0+t1);

/*	
    double cr = creal(c);
	double ci = cimag(c); 

	double func(double x) {
		return x/sqrt((x-a)*(x-b)*(x*x-2.*x*cr+cr*cr+ci*ci));
	}

	if (X1<a) error("integral_R_r2_cc: X1<a");	
	if (X1<b) error("integral_R_r2_cc: X1<b");	
	if (X1>X2) error("integral_R_r2_cc: X1>X2");	
		
	double qromb(double (*func)(double), double a, double b); 
	double x = qromb(&func, X1, X2);
	//if (x==0) fprintf(stderr,"a=%e  m=%e  j=%e  umax=%e\n",a,m,j,umax);
	return x;
*/
//??????????	
	
}


DEVICEFUNC
double integral_R_r2_re(double a, double b, double c, double d, double X)
// int_a^X x^2/sqrt((x-a)(x-b)(x-c)(x-d))
// X > a > b > c > d
// Eq. 258.11 (Byrd & Friedman)
{
	double m2 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt(((b-d)*(X-a))/((a-d)*(X-b))); 
	double u  = jacobi_isn(sn, m2);//elliptic_f_sin(sn,m);
	double a2 = (a-d)/(b-d);
	double b2 = ((a-d)*b)/(a*(b-d)); 
	double Z  = integral_Z2(a2,b2,u,m2) - integral_Z2(a2,b2,0,m2);
	return sqr(a)*2.0/sqrt((a-c)*(b-d)) * Z;
}


DEVICEFUNC
//double integral_R_r2_cc(double a, double b, double complex c, double X1, double X2)
double integral_R_r2_cc(double a, double b, sim5complex c, double X1, double X2)
// int_X1^X2 x^2/sqrt((x-a)(x-b)(x-c)(x-c*))
// X > a > b;  c = u + i*v;
// Eq. 260.03 (Byrd & Friedman)
{
	double u  = creal(c); 
	double v2 = sqr(cimag(c)); 
	double A  = sqrt(sqr(a-u) + v2); 
	double B  = sqrt(sqr(b-u) + v2); 
	double m  = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double g  = 1./sqrt(A*B);
	double alpha1 = (B*a+b*A)/(B*a-b*A);
	double alpha2 = (B+A)/(B-A);

	double u1 = elliptic_f_cos((X1*(A-B)+a*B-b*A)/(X1*(A+B)-a*B-b*A),m);
	double u2 = elliptic_f_cos((X2*(A-B)+a*B-b*A)/(X2*(A+B)-a*B-b*A),m);
	
	double t0 = pow(alpha1,2.)            * (integral_R0(u2,m)        - integral_R0(u1,m));
	double t1 = 2.*alpha1*(alpha2-alpha1) * (integral_R1(alpha2,u2,m) - integral_R1(alpha2,u1,m));
	double t2 = pow(alpha2-alpha1,2.)     * (integral_R2(alpha2,u2,m) - integral_R2(alpha2,u1,m));

	return pow((B*a-b*A)/(B+A),2.)*g * (t0+t1+t2);
/*
	double cr = creal(c);
	double ci = cimag(c); 

	double func(double x) {
		return x*x/sqrt((x-a)*(x-b)*(x*x-2.*x*cr+cr*cr+ci*ci));
	}

	if (X1<a) error("integral_R_r2_cc: X1<a");	
	if (X1<b) error("integral_R_r2_cc: X1<b");	
	if (X1>X2) error("integral_R_r2_cc: X1>X2");	
		
	double qromb(double (*func)(double), double a, double b); 
	double x = qromb(&func, X1, X2);
	//if (x==0) fprintf(stderr,"a=%e  m=%e  j=%e  umax=%e\n",a,m,j,umax);
	return x;
*/
	
//??????????	
	
}


DEVICEFUNC
double integral_R_rp_re(double a, double b, double c, double d, double p, double X)
// int_a^X 1/[(x-p)*sqrt((x-a)(x-b)(x-c)(x-d))]
// X > a > b > c > d
// Eq. 258.39 (Byrd & Friedman)
{
	double m2 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt(((b-d)*(X-a))/((a-d)*(X-b))); 
	double u1 = jacobi_isn(sn, m2);//elliptic_f_sin(sn,m2);
	double a2 = (a-d)/(b-d);
	double c2 = ((p-b)*(a-d))/((p-a)*(b-d)); 
	return -2.0/sqrt((a-c)*(b-d))/(p-a) * (integral_Z1(c2,a2,u1,m2) - integral_Z1(c2,a2,0,m2));
}


DEVICEFUNC
double integral_R_rp_re_inf(double a, double b, double c, double d, double p)
// int_a^infty 1/[(x-p)*sqrt((x-a)(x-b)(x-c)(x-d))]
// infinity > a > b > c > d
// Eq. 258.39 (Byrd & Friedman)
{
	double m2 = ((b-c)*(a-d))/((a-c)*(b-d));
	double sn = sqrt((b-d)/(a-d)); 
	double u1 = jacobi_isn(sn, m2);//elliptic_f_sin(sn,m2);
	double a2 = (a-d)/(b-d);
	double c2 = ((p-b)*(a-d))/((p-a)*(b-d)); 
	return -2.0/sqrt((a-c)*(b-d))/(p-a) * (integral_Z1(c2,a2,u1,m2) - integral_Z1(c2,a2,0,m2));
}


DEVICEFUNC
//double integral_R_rp_cc2(double a, double b, double complex c, double p, double X1, double X2)
double integral_R_rp_cc2(double a, double b, sim5complex c, double p, double X1, double X2)
// int_X1^X2 1/[(x-p)*sqrt((x-a)(x-b)(x-c)(x-d))]
// X1 > a > b; X1 > p; c = u + i*v;  d = c* 
// Eq. 260.04 (Byrd & Friedman)
{
    #ifndef CUDA
    if (X1<a) error("integral_R_rp_cc2: X1<a");	
    if (X1<b) error("integral_R_rp_cc2: X1<b");	
    if (X1<p) error("integral_R_rp_cc2: X1<p");
    #endif
    
	double u  = creal(c); 
	double v2 = sqr(cimag(c)); 
	double A  = sqrt(sqr(a-u) + v2); 
	double B  = sqrt(sqr(b-u) + v2); 
	double m  = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double g  = 1./sqrt(A*B);
	double alpha1 = (B*a+b*A-p*A-p*B)/(B*a-b*A+p*A-p*B);
	double alpha2 = (B+A)/(B-A);

//if (fabs(fabs(alpha1)-1.0)<1e-10) {fprintf(stderr,"WRN: alpha1==1 in integral_R_rp_cc2 (%.10e; a=%.8e p=%.8e )\n",alpha1,a,p); return 0.0;}

	double u1 = elliptic_f_cos((X1*(A-B)+a*B-b*A)/(X1*(A+B)-a*B-b*A),m);
	double u2 = elliptic_f_cos((X2*(A-B)+a*B-b*A)/(X2*(A+B)-a*B-b*A),m);

	double t0 = alpha2 * (integral_R0(u2,m) - integral_R0(u1,m));
	double t1 = (alpha1-alpha2) * (integral_R1(alpha1,u2,m) - integral_R1(alpha1,u1,m));

	return (B-A)*g/(B*a+b*A-p*A-p*B) * (t0+t1);
}


DEVICEFUNC
//double integral_R_rp_cc2_inf(double a, double b, double complex c, double p, double X1)
double integral_R_rp_cc2_inf(double a, double b, sim5complex c, double p, double X1)
// int_X1^inf 1/[(x-p)*sqrt((x-a)(x-b)(x-c)(x-d))]
// X1 > a > b; X1 > p; c = u + i*v;  d = c* 
// Eq. 260.04 (Byrd & Friedman)
{
    #ifndef CUDA
    if (X1<a) error("integral_R_rp_cc2_inf: X1<a");	
    if (X1<b) error("integral_R_rp_cc2_inf: X1<b");	
    if (X1<p) error("integral_R_rp_cc2_inf: X1<p");
    #endif
    
	double u  = creal(c); 
	double v2 = sqr(cimag(c)); 
	double A  = sqrt(sqr(a-u) + v2); 
	double B  = sqrt(sqr(b-u) + v2); 
	double m  = (sqr(A+B) - sqr(a-b))/(4.*A*B);
	double g  = 1./sqrt(A*B);
	double alpha1 = (B*a+b*A-p*A-p*B)/(B*a-b*A+p*A-p*B);
	double alpha2 = (B+A)/(B-A);

//if (fabs(fabs(alpha1)-1.0)<1e-10) {fprintf(stderr,"WRN: alpha1==1 in integral_R_rp_cc2_inf (%.10e; a=%.8e p=%.8e )\n",alpha1,a,p); return 0.0;}

	double u1 = elliptic_f_cos((X1*(A-B)+a*B-b*A)/(X1*(A+B)-a*B-b*A),m);
	double u2 = elliptic_f_cos((A-B)/(A+B),m);
	
	double t0 = alpha2 * (integral_R0(u2,m) - integral_R0(u1,m));
	double t1 = (alpha1-alpha2) * (integral_R1(alpha1,u2,m) - integral_R1(alpha1,u1,m));

	return (B-A)*g/(B*a+b*A-p*A-p*B) * (t0+t1);
}




//---------------------------------------------------------------------------

DEVICEFUNC INLINE
double integral_T_m0(double a2, double b2, double X)
// int_X^b 1/sqrt((a^2+x^2)(b^2-x^2))
// b > X >= 0 
// Eq. 213.00 (Byrd & Friedman)
{
	double m = b2/(a2+b2);
	return 1./sqrt(a2+b2) * jacobi_icn(X/sqrt(b2),m);//elliptic_f_cos(X/sqrt(b2),m); 
}


DEVICEFUNC INLINE
double integral_T_m2(double a2, double b2, double X)
// int_X^b x^2/sqrt((a^2+x^2)(b^2-x^2))
// b > X >= 0 
// Eq. 213.06 (Byrd & Friedman)
{
	double m  = b2/(a2+b2);
	double cn = X/sqrt(b2);
    return b2/sqrt(a2+b2) * (integral_C2_cos(cn,m) - integral_C2(0,m));
}


DEVICEFUNC INLINE
double integral_T_mp(double a2, double b2, double p, double X)
// int_X^b 1/[(p-x^2)*sqrt((a^2+x^2)(b^2-x^2))]
// -b <= X <= +b 
// Eq. 213.02 (Byrd & Friedman)
{
    #ifndef CUDA
	if (fabs(X) > sqrt(b2)) error("integral_T_mp: invalid input ((X<0)||(X>b))");
	if (p==b2) error("integral_T_mp: invalid input (p==b2)");//{fprintf(stderr,"p==b2\n"); return 0.0;}
    #endif
	double m = b2/(a2+b2);	
    double n = b2/(b2-p);

    if (X >= 0.0) 
        return 1./sqrt(a2+b2)/(p-b2) * elliptic_pi_cos(X/sqrt(b2), n, m);
    else
        return 1./sqrt(a2+b2)/(p-b2) * (2.*elliptic_pi_complete(n, m) - elliptic_pi_cos(-X/sqrt(b2), n, m));
}







