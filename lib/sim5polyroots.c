/*
#include "sim5config.h"
#ifndef CUDA
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
//#include <complex.h>
#include "sim5math.h"
#include "sim5polyroots.h"
#endif
*/

// Calculate the roots of z^2 + p z + q = 0, and the number of real roots. Here
// p and q are complex numbers. We use the algorithm recommended in "Numerical
// Recipes". Input: pr = Re(p), pi = Im (p), qr = Re(q), qi = Im(q). Output: 
// the number of real roots, and solutions of z. zr = Re(z), zi = Im(z).
DEVICEFUNC
int quadratic_eq(double pr, double pi, double qr, double qi, double *zr, double *zi)
{
/*
  complex bb, bb2, cc, del, qq, z1, z2;
  double tst;
  int i, nu_re=0;

  bb  = cmake(pr, +pi);
  bb2 = ccc(bb);
  cc  = cmake(qr, qi);

  del = csub(csqr(bb), cmulr(cc,4.0));
  tst = cre(cmulc(bb2,csqrt(del)));

  if (tst>=0.0)
    qq=cmulr(cadd(bb,csqrt(del)),-0.5);
  else
    qq=cmulr(csub(bb,csqrt(del)),-0.5);

  z1 = qq;
  z2 = cdiv(cc,qq);

  zr[0]=z1.re;
  zi[0]=z1.im;
  zr[1]=z2.re;
  zi[1]=z2.im;

  for (i=0;i<2;i++) {
	 if (zi[i]==0.0) ++nu_re;
  }

  return nu_re;
*/
  //double complex eye = I;
  //double complex bb, bb2, cc, del, qq, z1, z2;
  sim5complex bb, bb2, cc, del, qq, z1, z2;
  double tst;
  int i, nu_re=0;

  //bb=pr+eye*pi;
  //bb2=pr-eye*pi;
  //cc=qr+eye*qi;
  bb  = makeComplex(pr, pi);
  bb2 = makeComplex(pr, -pi);
  cc  = makeComplex(qr, qi);

  del=bb*bb-4.*cc;
  tst=creal(bb2*csqrt(del));

  if (tst>=0.0)
    qq=-0.5*(bb+csqrt(del));
  else
    qq=-0.5*(bb-csqrt(del));

  z1=qq;
  z2=cc/qq;

  zr[0]=creal(z1);
  zi[0]=cimag(z1);
  zr[1]=creal(z2);
  zi[1]=cimag(z2);

  for (i=0;i<2;i++)
    {
	 if (zi[i]==0.0)
	   ++nu_re;
    }
  return nu_re;
}


// Calculate the roots of z^3 + p z^2 + q z + r = 0, and the number of real roots. 
// Here p, q, and r are real numbers. We use the algorithm recommended in "Numerical
// Recipes". Input: p, q, r. Output: the number of real roots, and solutions of z. 
// zr = Re(z), zi = Im(z).
DEVICEFUNC
int cubic_eq(double p, double q, double r, double *zr, double *zi)
{
  double x1, x2, x3, y1, y2, y3;
  double theta, aa, bb, qq, rr;
  int sgnrr, i, nu_re=0;

  qq=(p*p-3.*q)/9.;
  rr=(2*p*p*p-9.*p*q+27*r)/54.;

  if (rr>=0.0)
    sgnrr=1;
  else
    sgnrr=-1;

  if ((rr*rr)<(qq*qq*qq))
    {
	 theta=acos(rr/sqrt(qq*qq*qq));

	 x1=-2*sqrt(qq)*cos(theta/3.)-p/3.;
	 x2=-2*sqrt(qq)*cos((theta+2*M_PI)/3.)-p/3.;
	 x3=-2*sqrt(qq)*cos((theta-2*M_PI)/3.)-p/3.;

	 y1=0.0;
	 y2=0.0;
	 y3=0.0;
    }
  else
    {
	 aa=-sgnrr*pow(fabs(rr)+sqrt(rr*rr-qq*qq*qq),1/3.);
	 if (aa!=0.0)
	   bb=qq/aa;
	 else
	   bb=0.0;

	 x1=(aa+bb)-p/3.;
	 x2=-0.5*(aa+bb)-p/3;
	 x3=x2;

	 y1=0;
	 y2=0.5*sqrt(3.)*(aa-bb);
	 y3=-y2;
    }

  zr[0]=x1;
  zi[0]=y1;
  zr[1]=x2;
  zi[1]=y2;
  zr[2]=x3;
  zi[2]=y3;

  for (i=0;i<3;i++)
    {
	 if (zi[i]==0.0)
	   ++nu_re;
    }
  return nu_re;
}


DEVICEFUNC
void sort_roots_re(double *r1, double *r2, double *r3, double *r4)
// sort to the order r1 > r2 > r3 > r4 
{
	double t;
	if (*r2 > *r1) { t=*r1; *r1=*r2; *r2=t; }
	if (*r3 > *r1) { t=*r1; *r1=*r3; *r3=t; }
	if (*r4 > *r1) { t=*r1; *r1=*r4; *r4=t; }
	if (*r3 > *r2) { t=*r2; *r2=*r3; *r3=t; }
	if (*r4 > *r2) { t=*r2; *r2=*r4; *r4=t; }
	if (*r4 > *r3) { t=*r3; *r3=*r4; *r4=t; }
}


// Sort real and complex numbers in the following way: 1) separate complex 
// numbers and real numbers, putting complex ones ahead of real ones; 2) sort 
// the real numbers in an descent order. r1=Re(r), r2=Im(r), s is the beginning
// index for real roots. E.g., for four numbers 2, 1+3i, 1-3i, -4, the subroutine
// sort_mix gives 1+3i, 1-3i, 2, -4 as output numbers, and s = 2. [Note, s 
// starts from 0.]
DEVICEFUNC
void sort_mix(double *r1, double *r2, int *s)
{
  int N=4;
  double rr1[N], rr2[N], tmp;
  int i, j, k=0;

  *s=0;

  for (i=0; i<N; i++) {
    if (r2[i] != 0.) {
        rr1[*s]=r1[i];
        rr2[*s]=r2[i];
        *s+=1;
	}
      else
	{
        rr1[N-1-k]=r1[i];
        rr2[N-1-k]=0.0;
        k+=1;
	}
  }

  for (i=0; i<N; i++){
      r1[i]=rr1[i];
      r2[i]=rr2[i];
  }

  for (i=*s; i<N; i++) {
      for (j=0; j<(N-i); j++) {
        if (r1[i+j]>r1[i]) {
            tmp=r1[i+j];
            r1[i+j]=r1[i];
            r1[i]=tmp;
	    }
      }
    }
}


// Sort real and complex numbers in the following way: 1) separate complex 
// numbers and real numbers, putting real ones ahead of complex ones; 2) sort 
// the real numbers in an descent order. r1=Re(r), r2=Im(r), s is the beginning
// index for complex roots. E.g., for four numbers 2, 1+3i, 1-3i, -4, the subroutine
// sort_mix gives 2, -4, 1+3i, 1-3i as output numbers, and s = 2. [Note, s 
// starts from 0.]
DEVICEFUNC
void sort_mix2(double *r1, double *r2, int *s)
{
  int N=4;
  double rr1[N], rr2[N], tmp;
  int i, j, k=0;

  *s=0;

  for (i=0; i<N; i++) {
    if (r2[i] == 0.) {
        rr1[*s]=r1[i];
        rr2[*s]=0.0;
        *s+=1;
	}
      else
	{
        rr1[N-1-k]=r1[i];
        rr2[N-1-k]=r2[i];
        k+=1;
	}
  }

  for (i=0; i<N; i++){
      r1[i]=rr1[i];
      r2[i]=rr2[i];
  }

  for (i=0; i<*s; i++) {
      for (j=0; j<(*s-i); j++) {
        if (r1[i+j]>r1[i]) {
            tmp=r1[i+j];
            r1[i+j]=r1[i];
            r1[i]=tmp;
	    }
      }
    }
}



DEVICEFUNC
//void sort_roots(int *s, double complex *z1, double complex *z2, double complex *z3, double complex *z4)
void sort_roots(int *s, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4)
// s returns no. of real roots
{
  int N=4;
  //double complex rr1[N], rr2[N], tmp;
  sim5complex rr1[N], rr2[N], tmp;
  int i, j, k=0;

  *s=0;
  
  rr1[0] = *z1;
  rr1[1] = *z2;
  rr1[2] = *z3;
  rr1[3] = *z4;

  for (i=0; i<N; i++) {
      if (cimag(rr1[i])==0.) {
          rr2[*s] = rr1[i];
          *s += 1;
      }
  }

  k=*s;
  for (i=0; i<N; i++) {
      if (cimag(rr1[i])!=0.) {
          rr2[k] = rr1[i];
          k += 1;
      }
  }

  for (i=0; i<*s; i++) {
      for (j=0; j<(*s-i); j++) {
        if (creal(rr2[i+j])>creal(rr2[i])) {
            tmp=rr2[i+j];
            rr2[i+j]=rr2[i];
            rr2[i]=tmp;
	    }
      }
    }
    
    *z1 = rr2[0];
    *z2 = rr2[1];
    *z3 = rr2[2];
    *z4 = rr2[3];
}



// Calculate the roots of z^4+ a3 z^3 + a2 z^2 + a1 z + a0 = 0, and the number of 
// real roots. Here a0, a1, a2, and a3 are real numbers. Input: a3, a2, a1. 
// Output: the number of real roots, and solutions of z. zr = Re(z), zi = Im(z).
DEVICEFUNC
int quartic_eq(double a3, double a2, double a1, double a0, double *zr, double *zi)
{
  double u1, pp, qq, rr, sup, del;
  double x1, x2, x3, x4, y1, y2, y3, y4;
  double zr2[2], zi2[2], zr3[3], zi3[3];
  int i, nu_re;

  pp=a2-3.*a3*a3/8.;
  qq=a1-0.5*a2*a3+0.125*pow(a3,3);
  rr=a0-0.25*a1*a3+a2*a3*a3/16.-3.*pow(a3,4)/256.;

  if (qq!=0)                     // then u1 \neq pp
    {
      cubic_eq(-pp,-4.*rr,4.*rr*pp-qq*qq,zr3,zi3);
      u1=zr3[0];                 // the first root given by cubic_eq must be real

      if (u1-pp>0.0)
	{
	  sup=sqrt(u1-pp);

	  quadratic_eq(sup,0,0.5*(u1-qq/sup),0,zr2,zi2);
	  x1=zr2[0];
	  x2=zr2[1];
	  y1=zi2[0];
	  y2=zi2[1];

	  quadratic_eq(-sup,0,0.5*(u1+qq/sup),0,zr2,zi2);
	  x3=zr2[0];
	  x4=zr2[1];
	  y3=zi2[0];
	  y4=zi2[1];
	}

      else
	{
	  sup=sqrt(pp-u1);

	  quadratic_eq(0,sup,0.5*u1,0.5*qq/sup,zr2,zi2);
	  x1=zr2[0];
	  x2=zr2[1];
	  y1=zi2[0];
	  y2=zi2[1];

	  quadratic_eq(0,-sup,0.5*u1,-0.5*qq/sup,zr2,zi2);
	  x3=zr2[0];
	  x4=zr2[1];
	  y3=zi2[0];
	  y4=zi2[1];
	}
    }
  else                           // if qq = 0
    {
      del=pp*pp-4.*rr;

      if (del>=0)
	{
	  quadratic_eq(0,0,0.5*(pp-sqrt(del)),0,zr2,zi2); 
	  x1=zr2[0];
	  x2=zr2[1];
	  y1=zi2[0];
	  y2=zi2[1];

	  quadratic_eq(0,0,0.5*(pp+sqrt(del)),0,zr2,zi2);
	  x3=zr2[0];
	  x4=zr2[1];
	  y3=zi2[0];
	  y4=zi2[1];
	}
      else
	{
	  quadratic_eq(0,0,0.5*pp,-0.5*sqrt(-del),zr2,zi2); 
	  x1=zr2[0];
	  x2=zr2[1];
	  y1=zi2[0];
	  y2=zi2[1];

	  quadratic_eq(0,0,0.5*pp,0.5*sqrt(-del),zr2,zi2);
	  x3=zr2[0];
	  x4=zr2[1];
	  y3=zi2[0];
	  y4=zi2[1];
	}
    }

  zr[0]=x1-0.25*a3;
  zi[0]=y1;
  zr[1]=x2-0.25*a3;
  zi[1]=y2;
  zr[2]=x3-0.25*a3;
  zi[2]=y3;
  zr[3]=x4-0.25*a3;
  zi[3]=y4;

  nu_re=0;
  for (i=0;i<4;i++)
    {
	 if (zi[i]==0.0)
	   ++nu_re;
    }

  if (nu_re == 4) sort_roots_re(&zr[0], &zr[1], &zr[2], &zr[3]);
  else
  if (nu_re == 2) sort_mix2(zr, zi, &i);

  return nu_re;
}


DEVICEFUNC
void quartic_eq_c(
	double a3, double a2, double a1, double a0, 
	int *nr, 
	//double complex *z1, double complex *z2, double complex *z3, double complex *z4)
	sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4)
{
	double zr[4];
	double zi[4];
	*nr = quartic_eq(a3,a2,a1,a0,zr,zi);
	*z1 = makeComplex(zr[0], zi[0]);
	*z2 = makeComplex(zr[1], zi[1]);
	*z3 = makeComplex(zr[2], zi[2]);
	*z4 = makeComplex(zr[3], zi[3]);
}




