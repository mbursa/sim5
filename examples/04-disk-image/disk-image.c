#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sim5lib.h"


const int image_dim_x = 1280;
const int image_dim_y =  720;



int main(int argc, char *argv[])
//***************************************************
{
    if (argc!=3) {
        fprintf(stderr, "Usage: %s <spin> <inclination>\n\n", argv[0]);
        fprintf(stderr, "spin - black hole rotation factor [0..0.999]\n");
        fprintf(stderr, "inclination - observer position [0..89]\n");
        exit(0);
    }

    // read paramters from command line
    double a   = (argc==3) ? atof(argv[1]) : 0.00;
    double inc = (argc==3) ? deg2rad(atof(argv[2])) : deg2rad(60.);
    
    if ((a<0.0) || (a>0.999) || (inc<0.0) || (inc>deg2rad(89.))) {
        fprintf(stderr, "ERROR: parameters out of range\n");
        exit(0);
    }

    fprintf(stderr, "Computing ... ");


    // allocate memory for images (2D arrays of float)
    // we have arrays for flux, g-factor
    float* image_f = (float*)calloc(image_dim_x*image_dim_y, sizeof(float));
    float* image_g = (float*)calloc(image_dim_x*image_dim_y, sizeof(float));

    // width of view - how big portion of the disk will the image cover
    double rms  = r_ms(a);
    double rmax = rms + 8.0;
    
    // setup disk model (Novikov-Thorne disk)
    disk_nt_setup(10.0, a, 0.1, 0.1, 0);
        
        
    // start the clock and calculate image
    clock_t cpu_t1, cpu_t2;
    cpu_t1 = clock();
    int ix, iy;
    // go over image and for each pixel independently to determine its brightness
    for (iy=0; iy<image_dim_y; iy++) for (ix=0; ix<image_dim_x; ix++) {
        // impact parameters
        // these are image coordinates scaled to rmax and shifted so 
        // that [0,0] is in the middle of the image
        double alpha = (((double)(ix)+.5)/(double)(image_dim_x)-0.5)*2.0*rmax;
        double beta  = (((double)(iy)+.5)/(double)(image_dim_y)-0.5)*2.0*rmax * ((double)image_dim_y/(double)image_dim_x);

        int error;
        double P,r,g,f;
        geodesic gd;

        // initialize photon trajectory for particular impact parameters [alpha,beta]
        geodesic_init_inf(inc, a, alpha, beta, &gd, &error);
        if (error) {
            fprintf(stderr, "photon rejected with error %d %e\n", error, gd.q);
            continue;
        }

        // find where the trajectory crosses equatorial plane and
        // obtain position parameter (returns NaN if there is no crossing)
        P = geodesic_find_midplane_crossing(&gd, 0);
        if (isnan(P)) continue;

        // from the position parameter get radius of disk intersection
        r = geodesic_position_rad(&gd, P);

        // if r<rms then the crossing lies inside the last stable orbit, 
        // where the disk cannot exist, thus radiation is zero.
        // outside of rms, determine disk radiation (f) and relativistic 
        // correction to it (g)
        if (r >= rms) {
            g = gfactorK(r, a, gd.l);
            f = disk_nt_flux(r);
            image_f[ix+image_dim_x*iy] = f*pow(g,4.);
            image_g[ix+image_dim_x*iy] = g;
            continue;
        }

        // the same thing for photons that make one orbit around the black hole
        // and may bring light from the bottom side of the disk
        P = geodesic_find_midplane_crossing(&gd, 1);
        if (isnan(P)) continue;

        r = geodesic_position_rad(&gd, P);

        if (r >= rms) {
            g = gfactorK(r, a, gd.l);
            f = disk_nt_flux(r);
            image_f[ix+image_dim_x*iy] = f*pow(g,4.);
            image_g[ix+image_dim_x*iy] = g;
            continue;
        }
    }
    // end if image calculation, read the clock
    cpu_t2 = clock();

    fprintf(stderr, "done\n");

    fprintf(stderr, "Profiling:\n");
    fprintf(stderr, "    photons: %d\n", image_dim_x*image_dim_y);
    fprintf(stderr, "    time: %.1f s\n", (double)(cpu_t2-cpu_t1)/CLOCKS_PER_SEC);
    fprintf(stderr, "    rate: %.1f photons/s\n", image_dim_x*image_dim_y/((double)(cpu_t2-cpu_t1)/CLOCKS_PER_SEC));

    // print image maps to stdout and free memory
    for (iy=0; iy<image_dim_y; iy++) {
        for (ix=0; ix<image_dim_x; ix++) printf("%d  %d  %e  %e\n", iy, ix, image_f[ix+image_dim_x*iy], image_g[ix+image_dim_x*iy]);
        printf("\n");
    }

    // free memory
    free(image_f);
    free(image_g);
    
    return 0;
}




