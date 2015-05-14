#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sim5lib.h"


#define EPS 1e-10


void test_raytrace();



int main() {

    test_raytrace();
    
    return 0;
}    




//---------------------------------------------------------------------------
// test methods 
//---------------------------------------------------------------------------


    
void test_raytrace() {
    srand (3000);
    
    int i;
    int N_total = 10000;
    int N_error = 0;
    int N_error_Q = 0;
    int N_error_KK = 0;
    int N_error_KF = 0;
    clock_t t1,t2;
    double time=0.0;
    double qq1,qq2,kk1,kk2,ff1,ff2,kf1,kf2;
    complex double wp1, wp2;
    
    FILE* log = fopen("test-raytrace.dat","w");   
    
    // raytrace N_total photons
    for (i=0; i<N_total; i++) {
        raytrace_data rtd;  // helper struct for raytrace
        sim5metric m;       // metric object
        sim5tetrad t;       // tetrad object
        int errors = 0;
    
        double bh_spin = rnd*0.999;
        double r_min = r_bh(bh_spin)*1.1;
        double r_max = 500.0;
    
        // set initial position
        double x[4];
        x[0] = 0.0;
        x[1] = r_min + 30.*rnd;
        x[2] = 2.*rnd-1.0;
        x[3] = 0.0;

        kerr_metric(bh_spin, x[1], x[2], &m);
        tetrad_zamo(&m, &t);

        // set initial direction (k-vector)
        double ang1 = rnd*M_2PI;
        double ang2 = rnd*M_PI;
        double n[4];
        double k[4];
        n[0] = -1.0;
        n[1] = sin(ang2)*cos(ang1);
        n[2] = sin(ang2)*sin(ang1);
        n[3] = cos(ang2);
        on2bl(n, k, &t);

        // initial polarization vector (f.k=0, f.f=1)
        double f[4], f_loc[4];
        double r1=rnd, r2=rnd, r3=rnd;
        f_loc[0] = 0.0;
        f_loc[1] = n[2]*r3 - n[3]*r2;
        f_loc[2] = n[3]*r1 - n[1]*r3;
        f_loc[3] = n[1]*r2 - n[2]*r1;
        on2bl(f_loc, f, &t);
        vector_norm_to(f, &m, 1.0);

        // checks of initial conditions
        kk1 = fabs(dotprod(k, k, &m));
        ff1 = fabs(dotprod(f, f, &m));
        kf1 = fabs(dotprod(f, k, &m));
        if (kk1>1e-10) fprintf(stderr, "ERR: KF>0 (%e)\n",kf1);
        if (kf1>1e-10) fprintf(stderr, "ERR: KF>0 (%e)\n",kf1);
        if (fabs(ff1-1.0)>1e-10) fprintf(stderr, "ERR: FF!=1 (%e)\n",ff1);

        // get motion constants
        qq1 = photon_carter_const(k, &m);
        wp1 = photon_wp_const(k, f, &m);

        // prepare raytrace
        raytrace_prepare(bh_spin, x, k, NULL, 0.01, RTOPT_NONE, &rtd);

        // do raytrace
        t1 = clock();
        while (1) {
            double dl = 1e9; // use maximal step
            raytrace(x, k, f, &dl, &rtd);
            // stop condition:
            if ((x[1] < r_min) || (x[1] > r_max)) break;
            // also stop if relative error this step is too large
            if (rtd.error>1e-3) break;
        }
        t2 = clock();

        // get total relative error        
        double error = raytrace_error(x, k, f, &rtd);
    
        time += (t2-t1)/(double)CLOCKS_PER_SEC;

        // construct polarization vector at a new position
        kerr_metric(bh_spin, x[1], x[2], &m);
        photon_polarization_vector(k, wp1, &m, f);

        // do checks
        wp2 = photon_wp_const(k, f, &m);
        qq2 = photon_carter_const(k, &m);
        kk2 = fabs(dotprod(k, k, &m));
        ff2 = fabs(dotprod(f, f, &m));
        kf2 = fabs(dotprod(k, f, &m));
        double QQ = fabs(qq2-qq1)/(qq1+1e-40);
        double FF = fabs(ff2-ff1)/(ff1+1e-40);
        double WP1 = fabs(creal(wp2)-creal(wp1))/(creal(wp1)+1e-40);
        double WP2 = fabs(cimag(wp2)-cimag(wp1))/(cimag(wp1)+1e-40);
        double WP  = fabs(cabs(wp2)-cabs(wp1))/(cabs(wp1)+1e-40);
        fprintf(log, "%d %e %e %e %e %e %e %e %e %e %e %e\n", rtd.pass, error, QQ, kk1, kk2, ff1, ff2, kf1, kf2, WP1, WP2, WP);
        if (error>1e-2) {errors|=1; N_error_Q++;}
        if (QQ  > 1e-03) {errors|=1; N_error_Q++;}
        //if (kk2 > 11e-03) {errors|=1; N_error_KK++;}
        if (kf2 > 1e-09) {errors|=1; N_error_KF++;}
        if (errors) N_error++;
    } //end of for

    printf("Total:  %d\n", N_total);
    printf("Time:   %.1f\n", time);
    printf("Photons/sec: %.1f\n", N_total/time);
    printf("Failed: %d\n", N_error);
    printf("Failed (KK): %d\n", N_error_KK);
    printf("Failed (KF): %d\n", N_error_KF);
    printf("Failed (Q):  %d\n", N_error_Q);
    
    fclose(log);
}



