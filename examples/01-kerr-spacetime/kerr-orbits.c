#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "sim5lib.h"


int main(int argc, char *argv[])
//***************************************************
{
    double a;
    double rbh, rph, rmb, rms;

    // print header
    fprintf(stdout, "# Locations of black-hole horizon (r_bh), photon orbit radius (r_ph),\n");
    fprintf(stdout, "# marginally bound orbit (r_mb) and marginally stable orbit (r_ms)\n");
    fprintf(stdout, "# in Kerr spacetime as a function of black-hole spin.\n");
    fprintf(stdout, "# Line format: spin  r_bh  r_ph  r_mb  r_ms\n");
    fprintf(stdout, "# ----\n");

    // go over spins from 0 to 1
    for (a=0.0; a<1.0; a+=0.01) {
        // get values of different radii
        rbh = r_bh(a);
        rph = r_ph(a);
        rmb = r_mb(a);
        rms = r_ms(a);

        // print it out
        fprintf(stdout, "%.4f  %e  %e  %e  %e\n", a, rbh, rph, rmb, rms);
    }

    return 0;
}

