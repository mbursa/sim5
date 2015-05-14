/*
#include "sim5config.h"
#ifndef CUDA
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim5config.h"
#include "sim5math.h"
#include "sim5kerr.h"
#include "sim5polyroots.h"
#include "sim5photontrace.h"
#endif
*/


#define frac_change(a,b) (fabs((a)/(b)-1.0))


static long raytrace_prof_N = 0;
static double raytrace_prof_t = 0.0;

void raytrace_profile() {
    fprintf(stderr, "raytrace_profile: N=%ld t=%.2e t/N=%.3e\n", 
        raytrace_prof_N, raytrace_prof_t, raytrace_prof_t/(double)(raytrace_prof_N));
}


DEVICEFUNC
void kerr_raytrace_prepare(double bh_spin, double x[4], double k[4], double dk[4], double f[4], raytrace_data* rtd)
{
    rtd->pass = 0;
    rtd->error = 0.0;

    sim5metric m;
    flat_metric(bh_spin, x[1], x[2], &m);
    rtd->kt = k[0]*m.g00 + k[3]*m.g03;
    rtd->kf = k[3]*m.g33 + k[0]*m.g03;
    rtd->Q  = photon_carter(&m, k);

    // evaluate 4-momentum derivative
    double G[4][4][4];
    kerr_connection(bh_spin, x[1], x[2], G);
    Gamma(G,k,dk);
    
    if (f) {
        rtd->K1 = 0.0;
        rtd->K2 = 0.0;
        // check f.k=0
        // evaluate K1,K2
    }
}



DEVICEFUNC
void kerr_raytrace(double bh_spin, double x[4], double k[4], double dk[4], double f[4], double *step, raytrace_data* rtd)
//***************************************************
// routine for geodesic integration: makes one step along the geodesic and 
// updates input vectors with new values
// inputs: <bh_spin> and vectors of position <x>, momentum <k>, momentum derivative <dk>,
//         polarization vector <f>, step size <dl> and <options>
// output: updates all vectors <x>, <k>, <dk>, <f>
// NOTE: polarization vector is optional, if it is NULL on input it is ignored
// see Dolence+2009 (http://adsabs.harvard.edu/abs/2009ApJS..184..387D)
{
    const float target_error = 1e-3;
    const float dl_reduction_factor = 5.0;

    int i;
    int iter = 0;
    double xp[4], kp[4], dkp[4], kp_prev[4];
    double kt=rtd->kt, kf=rtd->kf;
    double G[4][4][4];
    float k_frac_error;
    sim5metric m;
    
    // evaluate derivative of momentum from the geodesic equation; 
    // this function does the same thing as Gamma() from sim5kerr.c, except
    // in inline form it runs ~2 times faster (no need for passing G to the function)
    inline double k_deriv(int j, double _k[4]) {
        int a,b;
        double _dki = 0.0;
        for (a=0;a<4;a++) for (b=a;b<4;b++) _dki += -G[j][a][b]*_k[a]*_k[b];
        return _dki;
    }



    // limit step into a reasonable iterval
    double dl = minmax(1e-4, 5e0, *step) * dl_reduction_factor;
    
    rtd->pass++;
    
    do {
        iter++;
        if (iter>3) break;

        // reduce the time step by dl_reduction_factor in each pass 
        // for the first pass the reduction is taken into account in the initial value of dl
        dl /= dl_reduction_factor;
        
        // step 1: update of position (Dolence+09, Eq.14a)
        // we use m=cos(theta) for x[2], thus we need to add dm/d\theta=-sin(theta) factor
        xp[0] = x[0] + k[0]*dl + 0.5*dk[0]*sqr(dl);
        xp[1] = x[1] + k[1]*dl + 0.5*dk[1]*sqr(dl);
        xp[2] = x[2] +(k[2]*dl + 0.5*dk[2]*sqr(dl))*(-sqrt(1.-x[2]*x[2]));
        xp[3] = x[3] + k[3]*dl + 0.5*dk[3]*sqr(dl);

        // update metric
        kerr_metric(bh_spin, xp[1], xp[2], &m);

        //if ((xp[2]>+1.0) && (xp[2]<+1.0+1e-5)) xp[2] = +1.0-1e-5;
        //if ((xp[2]<-1.0) && (xp[2]>-1.0-1e-5)) xp[2] = -1.0+1e-5;
        if (fabs(xp[2])>1.0) {
            //fprintf(stderr,"#WRN: m>0 (%e dl=%.2e k=%.3e k=%.3e f=%.3e)\n",xp[2], dl, k[2]*dl,dk[2]*dl*dl, -sqrt(1.-xp[2]*xp[2]));
            //fprintf(stderr,"raytrace restart (%ld/%d): m>1\n",rtd->pass,iter);
            //dl /= 5.0;
            continue;
        }

        // step 2: update of momentum (Dolence+09, Eq.14b)
        for (i=0;i<4;i++) kp[i] = k[i] + dk[i]*dl;

        // presision check
        kt = kp[0]*m.g00 + kp[3]*m.g03;
        kf = kp[3]*m.g33 + kp[0]*m.g03;
        rtd->error = max(frac_change(kt,rtd->kt), frac_change(kf,rtd->kf));
        if (rtd->error > target_error*10.) {
            //fprintf(stderr,"raytrace restart (%ld/%d): step2\n",rtd->pass,iter);
            //dl /= 3.0;
            continue;
        }


        kerr_connection(bh_spin, xp[1], xp[2], G);

        // step 3: iteratively improve estimate of momentum based on its derivative
        int k_iter = 0;
        do {
            k_frac_error = 0.0;
            memcpy(kp_prev, kp, 4*sizeof(double)); //save value of kp to kp_prev
            
            for (i=0;i<4;i++) {
                // Dolence+09, Eq. (14c)
                dkp[i] = k_deriv(i,kp_prev);
        
                // Dolence+09, Eq. (14d)
                kp[i] = k[i] + 0.5*(dk[i] + dkp[i])*dl;

                k_frac_error += frac_change(kp[i],kp_prev[i]);
            }
        
            // quit, if either error in k is bellow treshold, or if k start to diverge
            if ((k_frac_error<target_error*1e-1) || (k_frac_error > 1.0)) break;

            k_iter++;
        } while(k_iter<3);

        // presision check: if k_frac_error is ~10 times target error, we better restart with lower step size
        // othervise we may try to continue with Runge-Kuta
        if (k_frac_error > target_error*1e+1) continue;
        
        if (k_frac_error > target_error*1e-1) {
            // if the above method has not achieved enough accuracy, 
            // we find the value of k^\mu using standard 4th-order Runge-Kuta method
            double kp1[4], kp2[4], kp3[4], kp4[4];
            for (i=0;i<4;i++) {
                kp1[i]  = k[i];
                kp2[i]  = k[i] + kp1[i]*dl/2.;
                kp3[i]  = k[i] + kp2[i]*dl/2.;
                kp4[i]  = k[i] + kp3[i]*dl;
            }
            double dkp1[4], dkp2[4], dkp3[4], dkp4[4];
            for (i=0;i<4;i++) {
                dkp1[i] = k_deriv(i,kp1);
                dkp2[i] = k_deriv(i,kp2);
                dkp3[i] = k_deriv(i,kp3);
                dkp4[i] = k_deriv(i,kp4);
            }
            for (i=0;i<4;i++) kp[i] = k[i] + dl/6.*(dkp1[i] + 2.*dkp2[i] + 2.*dkp3[i] + dkp4[i]);
        }

        // presision check
        kt = kp[0]*m.g00 + kp[3]*m.g03;
        kf = kp[3]*m.g33 + kp[0]*m.g03;
        rtd->error = max(frac_change(kt,rtd->kt), frac_change(kf,rtd->kf));


        /*
        //float Q1 = photon_carter(bh_spin, xp[1], xp[2], kp);
        //k_frac_change = max(k_frac_change, fabs(Q1-Q)/fabs(Q));
        //if (fabs(Q1-Q)/fabs(Q) > 1e-3) fprintf(stderr,"Carter Warning (%.4e)  dl=%.2e  Q=%.4e  Q1=%.4e\n", fabs(Q1-Q)/fabs(Q), dl, Q,Q1);
        double Q1 = photon_carter(bh_spin, xp[1], xp[2], k);
        if (fabs(Q1/rtd->Q-1.) > precision) {
            fprintf(stderr,"raytrace restart (%d r=%.3f)\n",iter, xp[1]);
            dl /= 5.0;
            iter++;
//            if (iter>10) fprintf(stderr,"raytrace iter too many (%d r=%.3f, prec=%.3e)\n",iter, xp[1], fabs(Q1-Q)/fabs(Q));
//            if (iter<=10) goto restart;
                continue;
        }
        */
    } while(rtd->error > target_error);


    // assign motion constants for k in the current step
    rtd->kt = kt;
    rtd->kf = kf;

    // assign x, k, and dk
    for (i=0;i<4;i++) {
        x[i]  = xp[i];
        k[i]  = kp[i];
        dk[i] = k_deriv(i,kp);
    }

    // assign the actual size of step taken
    *step = dl;
}




