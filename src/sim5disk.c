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
//! Loads an external library and binds its functions.


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dlfcn.h>



static void*  lib_handle = NULL;
static int    (*_diskmodel_init)(double M, double a, char* params) = NULL;
static void   (*_diskmodel_done)() = NULL;
static char*  (*_diskmodel_name)() = NULL;
static void   (*_diskmodel_help)() = NULL;
static double (*_diskmodel_r_min)() = NULL;
static double (*_diskmodel_flux)(double R) = NULL;
static double (*_diskmodel_lumi)() = NULL;
static double (*_diskmodel_mdot)() = NULL;
static double (*_diskmodel_temp)(double R) = NULL;
static double (*_diskmodel_sigma)(double R) = NULL;
static double (*_diskmodel_l)(double R) = NULL;
static double (*_diskmodel_vr)(double R) = NULL;
static double (*_diskmodel_h)(double R) = NULL;
static double (*_diskmodel_dhdr)(double R) = NULL;
static double (*_diskmodel_eval)(double R, int quantity) = NULL;
static void   (*_diskmodel_params)(FILE* output) = NULL;




int diskmodel_init(char *modellib, double M, double a, char* params)
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
    if (!bind_func((void**)&_diskmodel_help,  "diskmodel_help")) goto error;
    if (!bind_func((void**)&_diskmodel_params,"diskmodel_params")) goto error;
    if (!bind_func((void**)&_diskmodel_r_min, "diskmodel_r_min")) goto error;
    if (!bind_func((void**)&_diskmodel_mdot,  "diskmodel_mdot")) goto error;
    if (!bind_func((void**)&_diskmodel_lumi,  "diskmodel_lumi")) goto error;
    if (!bind_func((void**)&_diskmodel_flux,  "diskmodel_flux")) goto error;
    if (!bind_func((void**)&_diskmodel_temp,  "diskmodel_temp")) goto error;
    if (!bind_func((void**)&_diskmodel_sigma, "diskmodel_sigma")) goto error;
    if (!bind_func((void**)&_diskmodel_l,     "diskmodel_l")) goto error;
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
{
    _diskmodel_done();

    dlclose(lib_handle);
    lib_handle = NULL;
}




char* diskmodel_name()
{
    return _diskmodel_name();
}



double diskmodel_r_min()
{
    return _diskmodel_r_min();
}



double diskmodel_mdot()
{
    return _diskmodel_mdot();
}



double diskmodel_lumi()
{
    return _diskmodel_lumi();
}



double diskmodel_flux(double R)
{
    return _diskmodel_flux(R);
}


double diskmodel_temp(double R)
{
    return _diskmodel_temp(R);
}


double diskmodel_sigma(double R)
{
    return _diskmodel_sigma(R);
}


double diskmodel_l(double R)
{
    return _diskmodel_l(R);
}



double diskmodel_vr(double R)
{
    return _diskmodel_vr(R);
}



double diskmodel_h(double R)
{
    return _diskmodel_h(R);
}



double diskmodel_dhdr(double R)
{
    return _diskmodel_dhdr(R);
}


double diskmodel_eval(double R, int quantity)
{
    return _diskmodel_eval(R, quantity);
}



void diskmodel_params(FILE* output)
{
    _diskmodel_params(output);
}



void diskmodel_dump(char* filename)
{
    FILE* output = stderr;
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
            diskmodel_l(r),
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

