#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sim5lib.h"


#define EPS 1e-10


void test_raytrace();
void test_geodesic_init_src();
void test_ntdisk();
void test__gauss_distribution();
void test__interpolation();
void test__connection();


int main() {

    //test_ntdisk();

    //test__interpolation();
    //test__gauss_distribution();

    //test_raytrace();

    //test_geodesic_init_src();

    test__connection();


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

        double bh_spin = sim5urand()*0.999;
        double r_min = r_bh(bh_spin)*1.1;
        double r_max = 500.0;

        // set initial position
        double x[4];
        x[0] = 0.0;
        x[1] = r_min + 30.*sim5urand();
        x[2] = 2.*sim5urand()-1.0;
        x[3] = 0.0;

        kerr_metric(bh_spin, x[1], x[2], &m);
        tetrad_zamo(&m, &t);

        // set initial direction (k-vector)
        double ang1 = sim5urand()*PI;
        double ang2 = sim5urand()*PI2;
        double n[4];
        double k[4];
        n[0] = -1.0;
        n[1] = sin(ang2)*cos(ang1);
        n[2] = sin(ang2)*sin(ang1);
        n[3] = cos(ang2);
        on2bl(n, k, &t);

        // initial polarization vector (f.k=0, f.f=1)
        double f[4], f_loc[4];
        double r1=sim5urand(), r2=sim5urand(), r3=sim5urand();
        f_loc[0] = 0.0;
        f_loc[1] = n[2]*r3 - n[3]*r2;
        f_loc[2] = n[3]*r1 - n[1]*r3;
        f_loc[3] = n[1]*r2 - n[2]*r1;
        on2bl(f_loc, f, &t);
        vector_norm_to(f, 1.0, &m);

        // checks of initial conditions
        kk1 = fabs(dotprod(k, k, &m));
        ff1 = fabs(dotprod(f, f, &m));
        kf1 = fabs(dotprod(f, k, &m));
        if (kk1>1e-10) fprintf(stderr, "ERR: KF>0 (%e)\n",kf1);
        if (kf1>1e-10) fprintf(stderr, "ERR: KF>0 (%e)\n",kf1);
        if (fabs(ff1-1.0)>1e-10) fprintf(stderr, "ERR: FF!=1 (%e)\n",ff1);

        // get motion constants
        qq1 = photon_carter_const(k, &m);
        wp1 = polarization_constant(k, f, &m);

        // prepare raytrace
        raytrace_prepare(bh_spin, x, k, 0.01, RTOPT_NONE, &rtd);

        // do raytrace
        t1 = clock();
        while (1) {
            double dl = 1e9; // use maximal step
            raytrace(x, k, &dl, &rtd);
            // stop condition:
            if ((x[1] < r_min) || (x[1] > r_max)) break;
            // also stop if relative error this step is too large
            if (rtd.error>1e-3) break;
        }
        t2 = clock();

        // get total relative error
        double error = raytrace_error(x, k, &rtd);

        time += (t2-t1)/(double)CLOCKS_PER_SEC;

        // construct polarization vector at a new position
        kerr_metric(bh_spin, x[1], x[2], &m);
        polarization_vector(k, wp1, &m, f);

        // do checks
        wp2 = polarization_constant(k, f, &m);
        qq2 = photon_carter_const(k, &m);
        kk2 = fabs(dotprod(k, k, &m));
        ff2 = fabs(dotprod(f, f, &m));
        kf2 = fabs(dotprod(k, f, &m));
        double QQ = fabs(qq2-qq1)/(qq1+1e-40);
        //double FF = fabs(ff2-ff1)/(ff1+1e-40);
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



void test_geodesic_init_src()
{
    double a   = 0.099;
    double inc = deg2rad(170.);
    double rmax = 30.0;
    //double dr   = 0.1;

    int x, y;
    int N = 50;
    int tested = 0;

    for (y=0; y<N; y++) {
        for (x=0; x<N; x++) {
            //printf("------------\n");
            double alpha = (((double)(x)+.5)/(double)(N)-0.5)*2.0*rmax;
            double beta  = (((double)(y)+.5)/(double)(N)-0.5)*2.0*rmax;

            int pa, status;
            double P,r,phi1;
            geodesic gd1;
            geodesic gd2;

            geodesic_init_inf(inc, a, alpha, beta, &gd1, &status);
            if (!status) {
                //fprintf(stderr, "ERROR: cannot setup GD: %5d %5d  %e  %f\n", y, x, 0.0, 0.0);
                //printf("%.4f %.4f  %e  %e  %e\n", alpha, beta, 0.0, 0.0, -10000.);
                continue;
            }

            P = geodesic_find_midplane_crossing(&gd1, 0);
            if (isnan(P)) continue;
            pa = (P > gd1.Rpc);

            // from the position parameter get radius of disk intersection
            r = geodesic_position_rad(&gd1, P);
            if (isnan(r) || (r < r_bh(a))) continue;

            phi1 = geodesic_position_azm(&gd1, r, 0.0, P);

            //printf("%.4f %.4f  %e  %e\n", alpha, beta, r, phi);

            double k[4];
            photon_momentum(a, r, 0.0, gd1.l, gd1.q, pa?-1.0:+1.0, -1.0, k);


            geodesic_init_src(a, r, 0.0, k, pa, &gd2, &status);

            raytrace_data rtd;
            double x[4];
            vector_set(x, 0.0, r, 0.0, 0.0);
            raytrace_prepare(a, x, k, 0.01, RTOPT_NONE, &rtd);
            while (1) {
                double dl = 1e9; // use maximal step
                raytrace(x, k, &dl, &rtd);
                // stop condition:
                if ((x[1] < r_bh(a)) || (x[1] > 1e9)) break;
                // also stop if relative error this step is too large
                if (rtd.error>1e-2) break;
            }
            if (x[1] > 1e9) {
                //fprintf(stderr, "inc=%e/%e  phi=%e/%e\n", x[2], gd1.cos_i, reduce_angle_2pi(-x[3]), reduce_angle_2pi(phi1));
                fprintf(stderr, "d_inc=%+e  d_phi=%+e (%+.2f/%+.2f  %+.2e/%+.2e)\n", x[2]-gd1.cos_i, x[3]-phi1, acos(x[2])*180/M_PI, acos(gd1.cos_i)*180/M_PI, x[3], phi1);
                //x[2] vs cos_i
                //x[3] vs phi1;
            }


            tested++;
            if (fabs(gd1.cos_i-gd2.cos_i) > 1e-5) {
                printf("mp: %e / %e\n", gd1.m2p, gd2.m2p);
                printf("mm: %e / %e\n", gd1.m2m, gd2.m2m);
                printf("mK: %e / %e\n", gd1.mK, gd2.mK);
                printf("l: %e / %e\n", gd1.l, gd2.l);
                printf("q: %e / %e\n", gd1.q, gd2.q);
                printf("i: %e / %e\n", gd1.cos_i, gd2.cos_i);
            } else {
                //printf("ok\n");
            }
        }

        //printf("\n");
    }
    printf("tested: %d\n", tested);

}


void test_ntdisk()
{
    double mass  = 10.0;
    double spin  = 0.0;
    double mdot  = 1.0;
    double alpha = 0.1;
    disk_nt_setup(mass, spin, mdot, alpha, 0);
    disk_nt_lumi();
    printf("tested.\n");

}


void test__interpolation()
{
    const double x_min = -5.0;
    const double x_max = +5.0;
    const int N1 = 11;
    const int N2 = 100;
    
    int i;
    double x[N1];
    double y[N1];
    
    double gauss(double _x) { return exp(-sqr(_x)/2.)/sqrt(2*M_PI); }

    for (i=0; i<N1; i++) {
        x[i] = x_min + (double)(i)*(x_max-x_min)/(double)(N1-1);
        y[i] = gauss(x[i]);    
    }
    
    sim5interp in;
    sim5_interp_init(&in, x, y, N1, INTERP_DATA_COPY, INTERP_TYPE_SPLINE, 0);
    for (i=0; i<N2; i++) {
        double xx = x_min + (double)(i)*(x_max-x_min)/(double)(N2-1);
        double yy = sim5_interp_eval(&in, xx);
        printf("%.3e  %.3e  %.3e\n", xx, gauss(xx), yy);
    }
    sim5_interp_done(&in);
}



void test__gauss_distribution()
{
    const double x_min = -6.0;
    const double x_max = +6.0;
    const int N_pdf = 100;
    const int N_hits = 300000000;
    
    int i;
    sim5distrib d;
    
    double gauss_pdf(double _x) { return exp(-sqr(_x)/2.)/sqrt(2*M_PI); }

    distrib_init(&d, gauss_pdf, x_min, x_max, N_pdf);  
    printf("# norm=%e\n", d.norm);
    
    //for(i=0; i<d.icd.N; i++) printf("%e %e\n", d.icd.X[i], d.icd.Y[i]);
    //return;
 
    double* pdf_x = (double*)malloc(N_pdf*sizeof(double));
    int* pdf_y = (int*)malloc(N_pdf*sizeof(int));
    double dx = (x_max-x_min)/(double)(N_pdf-1);
    for (i=0; i<N_pdf; i++) {
        pdf_x[i] = x_min + (double)(i)*dx;
        pdf_y[i] = 0.0;
    }
    
    for (i=0; i<N_hits; i++) {
        //long ix = sim5_interp_search(pdf_x, distrib_hit(&d), 0, N_pdf-1);
        double x = distrib_hit(&d);
        int ix = (int)( ((x+dx/2.)-x_min)/dx );
        if ((ix>=0)&&(ix<N_pdf)) pdf_y[ix]++;
    }

    for (i=0; i<N_pdf; i++) printf("%.3e  %.3e  %.3e\n", pdf_x[i], (float)(pdf_y[i])/(float)(N_hits)/((x_max-x_min)/N_pdf), gauss_pdf(pdf_x[i])/d.norm);
    
    free(pdf_x);
    free(pdf_y);
    distrib_done(&d);
}



void test__connection()
{
    int N_total = 100000000;
    clock_t start_t, end_t;
    double time;
    double tmp = 1.0;

    void kerr_connection1(double a, double r, double m, double G[4][4][4]);

    start_t = clock();
    for (long i = 0; i<N_total; i++) {
        double G[4][4][4];
        double a = sim5urand()*0.999;
        double r = r_bh(a) + sim5urand()*20.0;
        double m = sim5urand();
        kerr_connection(a, r, m, G);
        tmp *= G[1][1][1];
    }
    end_t = clock();
    time = ((float)end_t - (float)start_t)/CLOCKS_PER_SEC;
}
