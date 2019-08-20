//************************************************************************
//    SIM5 library
//    sim5math.h - math macros
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************


//! \file sim5math.h
//! Mathematical macros.
//!
//! Some useful mathematical macros.

#ifndef _SIM5MATH_H
#define _SIM5MATH_H

#ifdef __cplusplus
extern "C" {
#endif


//#ifdef SIM5_SINGLE_PREC_MATH
//    #define double   float
//#else
//    #define double   double
//#endif


#define PI      3.14159265359                                   //!< PI
#define PI2     6.28318530718                                   //!< 2*PI
#define PI4     12.5663706144                                   //!< 4*PI
#define PI_half 1.57079632679                                   //!< PI/2

#define sqr(a)   ((a) * (a))                                    //!< quadratic power of a
#define sqr2(a)  ((a) * (a))                                    //!< quadratic power of a
#define sqr3(a)  ((a) * (a) * (a))                              //!< cubic power of a
#define sqr4(a)  ((a) * (a) * (a) * (a))                        //!< quartic power of a
#define sqrt3(a) cbrt(a)                                        //!< cubic root of a
#define sqrt4(a) pow(a,0.25)                                    //!< quartic root of a
#define odd(a) ((a%2==1)?1:0)                                   //!< is odd number
#define sign(a) ((a) >= 0.0 ? (+1.0) : (-1.0))                  //!< positive or negative
//#define hypot(a,b) (sqrt((a)*(a) + (b)*(b)))                    //!< 
#define deg2rad(a) ((a)/180.0*M_PI)                             //!< convert degrees to radians
#define rad2deg(a) ((a)*180.0/M_PI)                             //!< convert radians to degrees
#define EE(a) pow(10.0,a)                                       //!< power of 10
#define ave(a, b, w) ((1.0-(w))*(a) + (w)*(b))                  //!< weighted average of two values
#define logave(a, b, w) (exp((1.0-(w))*log(a) + (w)*log(b)))    //!< weighted logarithmic average of two values
#define inrange(a, min, max) (((a)>=(min))&&((a)<=(max)))       //!< value in range



#ifdef CUDA
    #define ComplexI    makeComplex(0.0,1.0)                    //!< complex unit (sqrt(-1))
    typedef double2     sim5complex;                            //!< proxy to built-in "double2" type
#else
    #include <complex.h>
    #undef I
    #define ComplexI _Complex_I                                 //!< complex unit (sqrt(-1))
    typedef double _Complex sim5complex;                         //!< proxy to built-in "complex double" type
#endif


//! \cond SKIP
DEVICEFUNC INLINE long sim5round(double num);

DEVICEFUNC INLINE long int factorial(long int n);

DEVICEFUNC INLINE double reduce_angle_pi(double phi);
DEVICEFUNC INLINE double reduce_angle_2pi(double phi);

DEVICEFUNC int ensure_range(double *val, double min, double max, double acc);


DEVICEFUNC void cartesian2spherical1(double x, double y, double z, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf);
DEVICEFUNC void cartesian2spherical2(double cos_h, double sin_f, double cos_f, double Vx, double Vy, double Vz, double* Vr, double* Vh, double* Vf);


DEVICEFUNC INLINE void sim5seed();
DEVICEFUNC INLINE unsigned long long sim5rand();
DEVICEFUNC INLINE double sim5urand();


DEVICEFUNC INLINE sim5complex makeComplex(double r, double i);
DEVICEFUNC INLINE sim5complex nullComplex();
//! \endcond

#ifdef __cplusplus
}
#endif


#endif
