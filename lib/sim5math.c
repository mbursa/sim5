//************************************************************************
//    SIM5 library
//    sim5math.h - math macros and routines
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************

/*
#include "sim5config.h"
#ifndef CUDA
#include <math.h>
#include "sim5math.h"
#endif
*/


DEVICEFUNC INLINE 
long sim5round(double num) {
 return (long)(num+0.5);
}


DEVICEFUNC INLINE 
long int factorial(long int n)
{
	if (n<=1) return(1);	else n=n*factorial(n-1);
	return(n);
}


DEVICEFUNC INLINE 
double reduce_angle_pi(double phi)
// returns angle in the interval [0..pi]
{
    while (phi < 0.0)  phi += 2.*M_PI;
    while (phi > M_PI) phi -= M_PI;
    return phi;
}


DEVICEFUNC INLINE 
double reduce_angle_2pi(double phi)
// returns angle in the interval [0..2pi]
{
    while (phi >= +M_2PI) phi -= M_2PI;
    while (phi <     0.0) phi += M_2PI;
    return phi;
}


DEVICEFUNC
int ensure_range(double *val, double min, double max, double acc)
{
	if (*val<min-acc) return 0;
	if (*val>max+acc) return 0;

	if (*val<min) *val = min;
	if (*val>max) *val = max;
	return 1;
}



#ifdef CUDA
    __device__ 
#endif 
long rndseed = 0x4d544750;
DEVICEFUNC 
double urand()
{
    #if LONG_MAX > (16807*2147483647)
    int const a    = 16807;      //ie 7**5
    int const m    = 2147483647; //ie 2**31-1
	rndseed = ((long)(rndseed * a))%m;
	return (double)(rndseed)/(double)(m);
    #else
    double const a    = 16807;      //ie 7**5
    double const m    = 2147483647; //ie 2**31-1
	double temp = rndseed * a;
	rndseed = (int) (temp - m * floor ( temp / m ));
	return (double)(rndseed)/m;
    #endif
}



void cartesian2spherical1(double x, double y, double z, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf)
// transforms cartesian vector [Vx,Vy,Vz] at point [x,y,z] to spherical coordinates [V_r,V_theta,V_phi]
// see: http://www.web-formulas.com/Math_Formulas/Linear_Algebra_Transform_from_Cartesian_to_Spherical_Coordinate.aspx
//      http://web.mit.edu/8.02t/www/materials/modules/ReviewB.pdf
{
    double r = sqrt(x*x + y*y + z*z);
    double cos_h = z/r;
    double sin_h = sqrt(1.-sqr(cos_h));
    double cos_f = x/r/sin_h;
    double sin_f = y/r/sin_h;
    (*Vr) = sin_h*cos_f*Vx + sin_h*sin_f*Vy + cos_h*Vz;
    (*Vh) = cos_h*cos_f*Vx + cos_h*sin_f*Vy - sin_h*Vz;
    (*Vf) =      -sin_f*Vx +       cos_f*Vy;
}

void cartesian2spherical2(double cos_h, double sin_f, double cos_f, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf)
// transforms cartesian vector [Vx,Vy,Vz] at point P=[r,theta,phi] to spherical coordinates [V_r,V_theta,V_phi]
// the point P is specified by cos_h=cos(theta), sin_f=sin)phi), cos_f=cos(phi)
// see: http://www.web-formulas.com/Math_Formulas/Linear_Algebra_Transform_from_Cartesian_to_Spherical_Coordinate.aspx
//      http://web.mit.edu/8.02t/www/materials/modules/ReviewB.pdf
{
    double sin_h = sqrt(1.-sqr(cos_h));
    (*Vr) = sin_h*cos_f*Vx + sin_h*sin_f*Vy + cos_h*Vz;
    (*Vh) = cos_h*cos_f*Vx + cos_h*sin_f*Vy - sin_h*Vz;
    (*Vf) =      -sin_f*Vx +       cos_f*Vy;
}




//---------------------------------------------------------------------
// complex algebra
//---------------------------------------------------------------------



DEVICEFUNC INLINE
sim5complex makeComplex(double r, double i)
{
#ifdef CUDA
    sim5complex res;
    res.x = r;
    res.y = i;
    return res;
#else
    sim5complex res;
    res = r + ComplexI*i;
    return res;
#endif
}


DEVICEFUNC INLINE
sim5complex nullComplex()
{
#ifdef CUDA
    return makeComplex(0.0,0.0);
#else
    return 0.0;
#endif
}


#ifdef CUDA

DEVICEFUNC INLINE
double creal (sim5complex a) 
{ 
    return a.x; 
}


DEVICEFUNC INLINE
double cimag (sim5complex a) 
{ 
    return a.y; 
}


DEVICEFUNC INLINE
sim5complex cconj (sim5complex a)
{
    return makeComplex(creal(a), -cimag(a));
}


DEVICEFUNC INLINE
double cabs(sim5complex a)
{
    // implementation guards against intermediate underflow and overflow by scaling 

    double x = creal(a);
    double y = cimag(a);
    double v, w, t;
    x = (double)fabs(x);
    y = (double)fabs(y);
    if (x > y) {
        v = x;
        w = y; 
    } else {
        v = y;
        w = x;
    }
    t = w / v;
    t = 1.0 + t * t;
    t = v * (double)sqrt(t);
    if ((v == 0.0) || (v > 3.402823466e38) || (w > 3.402823466e38)) {
        t = v + w;
    }
    return t;
}


DEVICEFUNC INLINE
sim5complex cadd(sim5complex a, sim5complex b)
{
    return makeComplex(creal(a)+creal(b), cimag(a)+cimag(b));
}


DEVICEFUNC INLINE
sim5complex csub(sim5complex a, sim5complex b)
{
    return makeComplex(creal(a)-creal(b), cimag(a)-cimag(b));
}


DEVICEFUNC INLINE
sim5complex cmul(sim5complex a, sim5complex b)
{
    sim5complex prod;
    prod = makeComplex(
        creal(a)*creal(b) - cimag(a)*cimag(b),
        creal(a)*cimag(b) + cimag(a)*creal(b)
    );
    return prod;
}


DEVICEFUNC INLINE
sim5complex cdiv(sim5complex a, sim5complex b)
{
    // implementation guards against intermediate underflow and overflow by scaling 

    sim5complex quot;
    double s = ((double)fabs((double)creal(b))) + 
               ((double)fabs((double)cimag(b)));
    double oos = 1.0 / s;
    double ars = creal(a) * oos;
    double ais = cimag(a) * oos;
    double brs = creal(b) * oos;
    double bis = cimag(b) * oos;
    s = (brs * brs) + (bis * bis);
    oos = 1.0 / s;
    quot = makeComplex(
        (ars*brs + ais*bis) * oos,
        (ais*brs - ars*bis) * oos
    );
    return quot;
}


DEVICEFUNC INLINE
sim5complex csqrt(sim5complex a)
{
	sim5complex c;
	double x,y,w,r;
	if ((a.x == 0.0) && (a.y == 0.0)) {
	    return nullComplex();
	} else {
		x=fabs(a.x);
		y=fabs(a.y);
		if (x >= y) {
			r=y/x;
			w=sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r=x/y;
			w=sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (a.x >= 0.0) {
			c.x=w;
			c.y=a.y/(2.0*w);
		} else {
			c.y=(a.y >= 0) ? w : -w;
			c.x=a.y/(2.0*c.y);
		}
		return c;
	}
}


DEVICEFUNC INLINE
sim5complex catan(sim5complex a)
{
/*
    double re = creal(a), im = cimag(a);
    sim5complex z;

    if (im == 0) return makeComplex(atan(re), 0.0);

    //This is a naive implementation which does not fully
    //take into account cancellation errors, overflow, underflow
    //etc.  It would benefit from the Hull et al treatment. 
    double r = hypot(re,im);
    double imag;
    double u = 2.*im / (1. + r*r);

    //the following cross-over should be optimized but 0.1 seems to work ok
    if (fabs(u) < 0.1) {
        imag = 0.25 * (log1p(u) - log1p(-u));
    }
    else {
        double A = hypot(re, im+1.);
        double B = hypot(re, im-1.);
        imag = 0.5*log(A/B);
    }

    if (re == 0.0) {
        if (I > 1) {
            z = makeComplex(M_PI_2, im);
        }
        else if (I < -1.0) {
            z = makeComplex(-M_PI_2, im);
        }
        else {
            z = makeComplex(0.0, im);
        };
    } else {
        z = makeComplex(0.5*atan2(2.*R, (1.+r)*(1.-r)), im);
    }

  return z;
*/
    return nullComplex();
}




 
inline DEVICEFUNC sim5complex operator+(sim5complex a, sim5complex b)
{
    return cadd(a,b);
}


inline DEVICEFUNC sim5complex operator+(sim5complex a, double b)
{
    return makeComplex(creal(a)+b, cimag(a));
}


inline DEVICEFUNC sim5complex operator+(double a, sim5complex b)
{
    return makeComplex(a+creal(b), cimag(b));
}


inline DEVICEFUNC void operator+=(sim5complex a, double b)
{
    a.x += b;
}


inline DEVICEFUNC void operator+=(sim5complex a, sim5complex b)
{
    a.x += b.x;
    a.y += b.y;
}


inline DEVICEFUNC sim5complex operator-(sim5complex a, sim5complex b)
{
    return csub(a,b);
}


inline DEVICEFUNC sim5complex operator-(sim5complex a, double b)
{
    return makeComplex(creal(a)-b, cimag(a));
}


inline DEVICEFUNC sim5complex operator-(double a, sim5complex b)
{
    return makeComplex(a-creal(b), -cimag(b));
}


inline DEVICEFUNC void operator-=(sim5complex a, double b)
{
    a.x -= b;
}


inline DEVICEFUNC void operator-=(sim5complex a, sim5complex b)
{
    a.x -= b.x;
    a.y -= b.y;
}

// negate
inline DEVICEFUNC sim5complex operator-(sim5complex &a)
{
    return makeComplex(-creal(a), -cimag(a));
}


inline DEVICEFUNC sim5complex operator*(sim5complex a, sim5complex b)
{
    return cmul(a,b);
}


inline DEVICEFUNC sim5complex operator*(double a, sim5complex b)
{
    return makeComplex(a*creal(b), a*cimag(b));
}


inline DEVICEFUNC sim5complex operator*(sim5complex a, double b)
{
    return makeComplex(b*creal(a), b*cimag(a));
}


inline DEVICEFUNC sim5complex operator/(sim5complex a, sim5complex b)
{
    return cdiv(a,b);
}


inline DEVICEFUNC sim5complex operator/(sim5complex a, double b)
{
    return makeComplex(creal(a)/b, cimag(a)/b);
}


#endif
