//************************************************************************
//    SIM5 library
//    sim5lib.h - library header file
//------------------------------------------------------------------------
//    Author:
//    Michal Bursa (bursa@astro.cas.cz)
//    Astronomical Institute
//    Bocni II 1401/1, 141-00 Praha 4, Czech Republic
//************************************************************************

#ifndef _SIM5LIB_H
#define _SIM5LIB_H

#include "sim5config.h"
#include "sim5include.h"

#ifndef CUDA
// external libs
#include "mt19937/mt19937.h"
#endif

// sim5lib parts
#include "sim5const.h"
#include "sim5math.h"
#include "sim5utils.h"
#include "sim5integration.h"

#ifndef CUDA
#include "sim5interpolation.h"
#include "sim5distributions.h"
#endif

#include "sim5roots.h"
#include "sim5elliptic.h"
#include "sim5polyroots.h"
#include "sim5raytrace.h"
#include "sim5kerr.h"
#include "sim5kerr-geod.h"

#ifndef CUDA
#include "sim5disk-nt.h"
#endif

#include "sim5polarization.h"
#include "sim5radiation.h"

#endif
