//************************************************************************
//    sim5integration.c - integration functions
//------------------------------------------------------------------------
//    Written by:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1a, 141-31 Praha 4, Czech Republic
//------------------------------------------------------------------------
//    Created: 09.12.2015
//************************************************************************


//! \file sim5integration.c
//! Numerical integration
//! 
//! Routines for simple numerical integrations of functions.


//! \cond SKIP
#define NMAX_TRAPEZOID 23
#define NMAX_SIMPSON 23
//! \endcond


//! \cond SKIP
DEVICEFUNC 
static void integrate_trapezoid_rule(double(*f)(double), double a, double b, int n, double *s)
//! Integration core routine based on trapezoid rule
//! - computes the `n`-th stage of refinement of an extended trapezoidal rule for function
//!   <f> and limits `a` and `b`
//! - when called with n=1, the routine returns `s` as the crudest estimate of \f$\int^a_b f(x) dx\f$
//! - subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy 
//!   of `s` by adding 2^(n-2) additional interior points.
//! - `s` must not be modified between sequential calls
//! - implementation is based on Numerical Recipes with the improvement by GSL,
//!   where a varible pointer is passed to the routine to store the value of 
//!   the latest refinement instead of that being a global variable
//! - calling scheme:
//!   `for(j=1; j<=M+1; j++) trapezoid_rule(func, a, b, j, &answer);`
{
    double x, tnm, sum, del;
    int it, j;

    if(n==1){
        *s = 0.5 * (b-a) * (f(b) + f(a));
    } else {
        for(it=1, j=1; j < n-1; j++)  it <<= 1;
        tnm = (double) it;
        del = (b-a) / tnm;
        x = a + 0.5 * del;
        for(sum=0.0, j=1; j<=it; j++, x+=del) { sum += f(x); }
        *s = 0.5 * (*s + del * sum);
    }
}
//! \endcond


DEVICEFUNC 
double integrate_trapezoid(double(*f)(double), double a, double b, double acc)
//! Integral of a function using trapezoid rule.
//! Computes the integral \f$ \int^a_b f(x) dx \f$ in a series of refinement steps until
//! the relative accuracy is better than requested or the maximum predefined number of steps is reached.
//!
//! Note: Accuracy should not be increased beyond ~10^-6 as roundoff errors start to accumulate 
//! if too many steps are taken.
//!
//! @param f pointer fo a function that is to be integrated
//! @param a interval of integration lower limit
//! @param b interval of integration upper limit
//! @param acc relative accuracy of integration
//!
//! @result Integral of the function over the interval [a,b].
{
    int n;
    double s = 0.0;
    double olds = DBL_MIN;  // any number that is unlikely to be the average of the function at its endpoints is ok

    for (n=1; n<=NMAX_TRAPEZOID; n++) {
        integrate_trapezoid_rule(f, a, b, n, &s);
        if (n > 3) {        // avoid spurious early convergence
            if ((fabs(s-olds) < acc*abs(olds)) || ((s==0.) && (olds==0.))) return s;
        }
        olds = s;
    }

    #ifndef CUDA
    warning("too many steps in integrate_trapezoid()");
    #endif

    return s;
}



DEVICEFUNC 
double integrate_simpson(double (*f)(double), double a, double b, double acc)
//! Integral of a function using Simpson rule.
//! Computes the integral \f$ \int^a_b f(x) dx \f$ in a series of refinement steps until
//! the relative accuracy is better than requested or the maximum predefined number of steps is reached.
//! Simpson rule is generally more efficient than trapezoid rule when the function to be integrated
//! has a finite 4th derivative (continuous 3rd  derivative).
//!
//! Note: Accuracy should not be increased beyond ~10^-6 as roundoff errors start to accumulate 
//! if too many steps are taken.
//!
//! @param f pointer fo a function that is to be integrated
//! @param a interval of integration lower limit
//! @param b interval of integration upper limit
//! @param acc relative accuracy of integration
//!
//! @result Integral of the function over the interval [a,b].
{
    int n;
    double s, st=0.0, ost, os;

    ost = os = -1.e50;

    for (n=1; n<=NMAX_SIMPSON; n++){
        integrate_trapezoid_rule(f, a, b, n, &st);
        s = (4.*st - ost) / 3.;
        if (n > 3) {        // avoid spurious early convergence
            if ((fabs(s-os) < acc*fabs(os)) || ((s==0.) && (os==0.))) return s;
        }
        os = s;
        ost = st;
    }
  
    #ifndef CUDA
    warning("too many steps in integrate_simpson()");
    #endif

    return s;
}


//-------------------------------------------------------------------------------------

//! \cond SKIP
double gammln(double xx)
// natural log of the complete Gamma function.
// returns the log of the gamma function for X
{
  double x,tmp,ser;
  static double cof[6]={76.18009173,-86.50532033,24.01409822,
			-1.231739516,0.120858003e-2,-0.536382e-5};
  int j;
  
  x=xx-1.0;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.0;
  for (j=0;j<=5;j++) {
    x += 1.0;
    ser += cof[j]/x;
  }
  return -tmp+log(2.50662827465*ser);
}
//! \endcond



DEVICEFUNC 
void gauleg(double x1, double x2, double x[], double w[], int n)
//! Gauss-Legendre quadrature.
//! Given the lower and upper limits 
//! of integration \f$x_1\f$ and \f$x_2\f$, and given \f$n\f$ points, this routine returns arrays 
//! \f$x[n]\f$ and \f$w[n]\f$ of length \f$n\f$, containing the abscissas and weights of 
//! the Gauss-Legendre n-point quadrature formula.
//!
//! To integrate a function \f$f(x)\f$ over the interval \f$[x_1,x_2]\f$, just do a simple summation:
//!
//!     gauleg(x1, x2, x, w, n)
//!     integral = sum(i=0, i<n-1, f(x[i])*w[i])
//!
//! Using n points, the method exactly integrates polynomials up to (2*n-1)'th degree.
//! If the function \f$f(x)\f$ is well approximated by a polynomial, the integral
//! will be very accurate.
//!
//! @param x1 interval of integration lower limit
//! @param x2 interval of integration upper limit
//! @param x pointer to an array where abscissas should be written
//! @param w pointer to an array where waights should be written
//! @param n number of quandrature points and also the minimal dimensions of x[] and w[] arrays
//!
//! @result Abscissas and weights are returned in x[] and w[] arrays.
{
    const double EPS = 3.0e-11;
	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1;

	m=(n+1)/2;
	xm=0.5*(x2+x1);
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++) {
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++) {
				p3=p2;
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
			}
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (fabs(z-z1) > EPS);
		x[i-1]=xm-xl*z;
		x[n+1-i-1]=xm+xl*z;
		w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
		w[n+1-i-1]=w[i-1];
	}
}



float qgaus(float (*func)(float), float a, float b)
//! Gauss-Legendre integral of a function.
//! Returns the integral of the function over the interval [a,b] computed using the by ten-point 
//! Gauss-Legendre integration: the function is evaluated exactly ten times at interior points 
//! in the range of integration.
//!
//! @param func pointer to a function to be integrated
//! @param a interval of integration lower limit
//! @param b interval of integration upper limit
//!
//! @result Value of the integral over the interval [a,b].
{
    int j;
    float xr,xm,dx,s;
    static float x[]={0.0,0.1488743389,0.4333953941,0.6794095682,0.8650633666,0.9739065285}; //The abscissas and weights First value of each array not used
    static float w[]={0.0,0.2955242247,0.2692667193,0.2190863625,0.1494513491,0.0666713443};
    xm=0.5*(b+a);
    xr=0.5*(b-a);
    s=0;  // Will be twice the average value of the function, since the
    for (j=1;j<=5;j++) {//    ten weights (five numbers above each used twice)
        dx=xr*x[j];  //  sum to 2.
        s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
    }
    return s *= xr;  //Scale the answer to the range of integration.
}



double qgaus_general(double w[], double y[], int N, double a, double b)
//! A generalized version of `qgauss()` function.
//! It gives the integral of the function over the interval [a,b] on a supplied 
//! grid of functional values y and ther weights.
//!
//! @param func pointer to a function to be integrated
//! @param w array of weights
//! @param y array of functional values y=f(x[i])
//! @param N length of w[] and y[] arrays
//! @param a interval of integration lower limit
//! @param b interval of integration upper limit
//!
//! @result Value of the integral over the interval [a,b].
{
    int j;
    double s  = 0.0;
    for (j=0; j<N; j++) s += w[j]*y[j];
    return s*(b-a);  //Scale the answer to the range of integration.
}


#undef NMAX_TRAPEZOID
#undef NMAX_SIMPSON

