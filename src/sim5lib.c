//************************************************************************
//    SIM5 library
//    sim5lib.c - libraray implementation file
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************

#include "sim5lib.h"

#ifndef CUDA
// external libs
#include "mt19937/mt19937.c"
#endif


// sim5lib parts
#include "sim5math.c"
#include "sim5utils.c"
#include "sim5integration.c"

#ifndef CUDA
#include "sim5interpolation.c"
#include "sim5distributions.c"
#endif

#include "sim5roots.c"
#include "sim5elliptic.c"
#include "sim5polyroots.c"
#include "sim5raytrace.c"
#include "sim5kerr.c"
#include "sim5kerr-geod.c"

#ifndef CUDA
#include "sim5disk-nt.c"
#endif

#include "sim5radiation.c"


