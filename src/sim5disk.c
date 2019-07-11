//************************************************************************
//    sim5disk.c
//------------------------------------------------------------------------
//    Date   : 2.10.2014
//    Author : Michal Bursa
//    e-mail : bursa@astro.cas.cz
//------------------------------------------------------------------------
//    (C) 2017 Michal Bursa
//************************************************************************

//! \file sim5disk.c
//! Wrapper for dynamic linking an external library with a disk model
//! 
//! Loads an external library and binds its functions. An external library (in Linux) is a compiled 
//! shared library with .so extension.
//!
//! The external library has to declare at least the bellow listed methods and provide an implementation 
//! of them. Calls to those methods in SIM5 are passed to the linked library and the result is returned 
//! by a wrapper function. See a demo of how to use this functionality in the examples folder.
//! 
//! NOTE: This unit uses static variables to store persistent information about the linked library. 
//! As a result, routines in this module are NOT thread-safe in a sense different threads cannot each 
//! link a different library. They can, however, all make calls to the already linked library.
//! For the same reasons, the routines declared here are not available to CUDA.



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>



static void*  lib_handle = NULL;
static int    (*_diskmodel_init)(double M, double a, char* params) = NULL;
static void   (*_diskmodel_done)() = NULL;
static char*  (*_diskmodel_name)() = NULL;
//static void   (*_diskmodel_help)() = NULL;
static double (*_diskmodel_r_min)() = NULL;
static double (*_diskmodel_flux)(double R) = NULL;
static double (*_diskmodel_lumi)() = NULL;
static double (*_diskmodel_mdot)() = NULL;
static double (*_diskmodel_sigma)(double R) = NULL;
static double (*_diskmodel_ell)(double R) = NULL;
static double (*_diskmodel_vr)(double R) = NULL;
static double (*_diskmodel_h)(double R) = NULL;
static double (*_diskmodel_dhdr)(double R) = NULL;
static double (*_diskmodel_eval)(double R, int quantity) = NULL;
static void   (*_diskmodel_params)(FILE* output) = NULL;




int diskmodel_init(char *modellib, double M, double a, char* params)
//! External disk model initialization.
//! Loads the compiled library, links is into the program and calls its initialization routine. 
//! All the required functions have to be declared in the library.
//!
//! @param modellib  filesystem path to the library (string)
//! @param M       mass of BH [Msun]
//! @param a       spin of BH [0..1]
//! @param params  parametres that are passed to the library initialization function
//! 
//! @result Return the result of the library's initialization function or -1 of the library could not 
//! loaded or be initialized.
{
    if (lib_handle) {
        fprintf(stderr, "diskmodel_init: cannot open '%s', close current model first\n", modellib);
        return -1;
    }

    lib_handle = dlopen(modellib, RTLD_NOW);
    if (!lib_handle) {
        fprintf(stderr, "diskmodel_init: opening failed (%s)\n", dlerror());
        return -1;
    }
    
    int bind_func(void** func_var, char* func_name) {
        char *errstr;
        (*func_var) = dlsym(lib_handle, func_name);
        errstr = dlerror();
        if (errstr) fprintf(stderr, "diskmodel_init: error in function binding (%s; %s)\n", func_name, errstr);
        if (!*func_var) fprintf(stderr, "diskmodel_init: missing function (%s)\n", func_name);
        return *func_var != NULL;
    }

    if (!bind_func((void**)&_diskmodel_init,  "diskmodel_init")) goto error;
    if (!bind_func((void**)&_diskmodel_done,  "diskmodel_done")) goto error;
    if (!bind_func((void**)&_diskmodel_name,  "diskmodel_name")) goto error;
    //if (!bind_func((void**)&_diskmodel_help,  "diskmodel_help")) goto error;
    if (!bind_func((void**)&_diskmodel_params,"diskmodel_params")) goto error;
    if (!bind_func((void**)&_diskmodel_r_min, "diskmodel_r_min")) goto error;
    if (!bind_func((void**)&_diskmodel_mdot,  "diskmodel_mdot")) goto error;
    if (!bind_func((void**)&_diskmodel_lumi,  "diskmodel_lumi")) goto error;
    if (!bind_func((void**)&_diskmodel_flux,  "diskmodel_flux")) goto error;
    if (!bind_func((void**)&_diskmodel_sigma, "diskmodel_sigma")) goto error;
    if (!bind_func((void**)&_diskmodel_ell,     "diskmodel_ell")) goto error;
    if (!bind_func((void**)&_diskmodel_vr,    "diskmodel_vr")) goto error;
    if (!bind_func((void**)&_diskmodel_h,     "diskmodel_h")) goto error;
    if (!bind_func((void**)&_diskmodel_dhdr,  "diskmodel_dhdr")) goto error;
    if (!bind_func((void**)&_diskmodel_eval,  "diskmodel_eval")) goto error;

    return _diskmodel_init(M, a, params);

    // error handling part    
    error:
    dlclose(lib_handle);
    lib_handle = NULL;
    return -1;
}




void diskmodel_done()
//! External disk model finitialization.
//! Frees memory and unlinks the libraray.
{
    _diskmodel_done();

    dlclose(lib_handle);
    lib_handle = NULL;
}




char* diskmodel_name()
//! Model name.
//! Returns a pointer to a string with the model's name.
{
    return _diskmodel_name();
}



double diskmodel_r_min()
//! Minimal radius of the disk (disk inner edge).
//! Gives minimal value for radius for which the functions provide valid results. 
//! E.g. for NT disk, this corresponds to the radius of the marginally stable orbit.
//!
//! @result Radius of disk inner edge [GM/c2]
{
    return _diskmodel_r_min();
}



double diskmodel_mdot()
//! Mass accretion rate.
//! Returns mass accretion rate in Eddington units of (Mdot_Edd*M).
//!
//! @result Mass accretion rate in Eddington units.
{
    return _diskmodel_mdot();
}



double diskmodel_lumi()
//! Total disk luminosity.
//! Luminosity is obtained by integrating local flux over the surface area of the disk (both sides)
//! going into the whole sky (4pi solid angle). 
//! The integration makes a proper transformation of the flux from local to coordinate frame, but 
//! it ignores other relativistic effects, e.g. light bending.
//!
//! \f[L = 2 * 2\pi \int F(r) (-U_t) r dr\f]
//!
//! @result Total disk luminosity of both surfaces [erg s-1]
{
    const float disk_rmax = 1e5;

    // integrate disk luminosity from r_ms to disk_rmax rg
    // - the integration uses 'logarithmic rule': L = \int f(x) dx \int f(x)*x d(log(x))

    double func_luminosity(double log_r)
    {
        double r = exp(log_r);
        // calculate U_t
        double gtt = -1. + 2./r;
        double gtf = -2.*disk_nt_bh_spin/r;
        double gff = sqr(r) + sqr(disk_nt_bh_spin) + 2.*sqr(disk_nt_bh_spin)/r;
        double Omega = 1./(disk_nt_bh_spin + pow(r,1.5));
        double U_t = sqrt(-1.0/(gtt + 2.*Omega*gtf + sqr(Omega)*gff)) * (gtt + Omega*gtf);
        double F = _diskmodel_flux(r);
        // dL = 2pi*r*F(r) dr, extra r comes from log integration
        return 2.*M_PI*r*2.0*(-U_t)*F * r;
    }

    double L = integrate_simpson(func_luminosity, log(disk_nt_disk_rms), log(disk_rmax), 1e-5);

    // fix units to erg/s
    L *= sqr(disk_nt_bh_mass*grav_radius);

    return L/(L_Edd*disk_nt_bh_mass);
}



double diskmodel_flux(double R)
//! Local flux from one side of the disk.
//! Provides radial radiation flux dependence measured in local frame, i.e. flux measured by an 
//! observer that is at rest with respect to the fluid.
//!
//! @param R radius of emission [GM/c2]
//!
//! @result Total outgoing flux from unit area on one side of the disk [erg cm-2 s-1]. 
{
    return _diskmodel_flux(R);
}


double diskmodel_sigma(double R)
//! Column density.
//! Returns midplane column density of the fluid, i.e. the fluid density integrated from midplane to the disk surface, 
//! at a given radius.
//!
//! @param R radius (measured in equatorial plane) [rg]
//!
//! @result Midplane column density in [g/cm2].
{
    return _diskmodel_sigma(R);
}


double diskmodel_ell(double R)
//! Specific angular momentum.
//! Returns specific angular momentum of the fluid at given radius. 
//!
//! @param R radius (measured in equatorial plane) [rg]
//!
//! @result Specific angular momentum in [g.u.].
{
    return _diskmodel_ell(R);
}



double diskmodel_vr(double R)
//! Radial velocity.
//! Returns bulk radial velocity of the fluid at given radius as measured by aan observer in the co-rotating frame.
//!
//! @param R radius (measured in equatorial plane) [rg]
//!
//! @result Radial velocity in [speed_of_light].
{
    return _diskmodel_vr(R);
}



double diskmodel_h(double R)
//! Surface height.
//! Returns the scale-height of the surface of the disk above midplane at given radius. 
//!
//! @param R radius (measured in equatorial plane) [rg]
//!
//! @result Scale-height [rg].
{
    return _diskmodel_h(R);
}



double diskmodel_dhdr(double R)
//! Derivative of surface height.
//! Returns surface profile as derivative \f$dH/dR\f$ of its height above midplane at given radius. 
//!
//! @param R radius (measured in equatorial plane) [rg]
//!
//! @result Derivative of surface height.
{
    return _diskmodel_dhdr(R);
}


double diskmodel_eval(double R, int quantity)
//! Other quantity evaluation.
//! Returns the value of a given quantity. The disk model may provide more quantities than the standard set. 
//! Additional quantities may be accessed using this function by providing the quantity identifier.
//!
//! @param R radius (measured in equatorial plane) [rg]
//! @param quantity quantity identified code
//!
//! @result The value of the requested quantity.
{
    return _diskmodel_eval(R, quantity);
}



void diskmodel_params(FILE* stream)
//! Prints model parameters.
//! Writes down the parameters of the model to the given stream.
 /!
//! @param stream stream to write to

    if (stream) _diskmodel_params(stream);
}



void diskmodel_dump(char* filename)
//! Prints the disk structure as a function of radius.
//! The function prints the profile of all quantities as a function of radius 
//! from r_ms to some outer radius (~2000 rg). It prints to a file identified by its path (overwrites existing) and 
//! if that is empty it prints to STDOUT.
//!
//! @param filename Path to a file that should be written. If NULL then it prints to STDOUT. 
{
    FILE* output = stdout;
    if (filename) output = fopen(filename, "w");
    if (!output) {
        fprintf(stderr, "diskmodel_dump: cannot open output (%s)\n", filename); 
        return;
    }

    fprintf(output, "# Disk model dump\n"); 
    fprintf(output, "#-------------------------------------------\n"); 
    fprintf(output, "# parameters:\n"); 
    diskmodel_params(output);
    fprintf(output, "#-------------------------------------------\n"); 
    fprintf(output, "# col1: radius [GM/c2]\n");
    fprintf(output, "# col2: flux (one side) [erg s-1 cm-2]\n");
    fprintf(output, "# col3: sigma [g cm-2]\n");
    fprintf(output, "# col4: specific angular momentum [none]\n");
    fprintf(output, "# col5: radial velocity [speed of light]\n");
    fprintf(output, "# col6: disk height [GM/c2]\n");
    fprintf(output, "# col7: disk slope (derivative dH/dR of height with respect to equatorial radius) [none]\n");
    fprintf(output, "#-------------------------------------------\n"); 
    double r;
    double r_min = diskmodel_r_min();
    for (r=r_min; r<5000.; r*=1.05) {
        printf("%e  %e  %e  %e  %+e  %+e  %+e\n",
            r, 
            diskmodel_flux(r),
            diskmodel_sigma(r),
            diskmodel_ell(r),
            diskmodel_vr(r),
            diskmodel_h(r),
            diskmodel_dhdr(r)
        );
    }
    fflush(output);
    if (filename) fclose(output);
}


/*
int main() {
    
    diskmodel_init("ntdisk/diskmodel-nt.so", 10.0, 0.0, NULL);
    if (!lib_handle) return -1;
    diskmodel_dump(NULL);
    diskmodel_done();

    return 0;
}
*/

