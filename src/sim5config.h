//************************************************************************
//    SIM5 library
//    sim5config.h - compile-time configuration options for SIM5
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************


#ifndef _SIM5CONFIG_H
#define _SIM5CONFIG_H



#ifdef __CUDACC__
    #define CUDA

    #define DEVICEFUNC   __device__
    #define HOSTFUNC     __host__
    #define INLINE       __inline__
#else
    #define DEVICEFUNC
    #define HOSTFUNC
    #define INLINE       inline   // extern in new in GNU99 standard (in GCC 5+)
#endif


#endif
