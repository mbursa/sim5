
# SIM5 - a general purpose library for GR raytracing and radiation transport

## Introduction

SIM5 is a C library with a Python interface that contains a collection of routines for relativistic raytracing and radiation transfer in GR. It has a special focus on raytracing from accretion disks, tori, hot spots or any other 3D configuration of matter in Kerr geometry, but it can be used with any other metric as well. It can handle both optically thick and thin sources as well as transport of polarization of the radiation and helps to calculate the propagation of light rays from the source to an observer through a curved spacetime. It supports parallelization and runs on GPUs.

## Features

The library contains
* a selection of physical and unit conversion constants related to black hole disk accretion
* an implementation of relativistic thin disk model
* basic routines related to GR and Kerr geometry in particular
* routines for computing null geodesics and GR raytracing with polarization 

The library is thread-safe and supports parallelization of the calculations at both CPUs (OpenMP/MPI) and GPUs (CUDA/OpenCL). 

## Installation

Acquire the source code by cloning the git repository

    git clone https://github.com/mbursa/sim5

The C code does not have any external library dependencies except the standard system ones. Typically, it is enough to run

    make

The compilation process creates a `lib` folder and puts two files there, `sim5lib.c`, `sim5lib.h` and `sim5lib.o`, which contain a complete C source code merged by the compiler process from all the individual source files in `src` folder. Use those files to statically link SIM5 in your code (most typical usage in C/C++ codes), see an example bellow.

### Compilation of the Python interface

The Python interface does have a few external dependecies and thus it does not compile by default when calling `make`. Before compiling the Python interface, make sure you have the following pre-requisities installed:

```bash
    # python headers
    apt install python-dev
    # SWIG - a tool that connects C/C++ programs with high-level programming languages
    apt install swig
```

Having those do

```bash
    make python
```

<!--
The third file is a compiled library used for run-time linking. It is used by the Python interface, but you may use it for dynamical linking of SIM5 to your code too.
-->

## Usage

The library can be directly used in C/C++ projects and with a proper interface it can be used also in other languages (Fortran, ...). It also works as a Python package making most of its functions accessible from Python scripts with no extra effort. 

### Using SIM5 in your C/C++ project (static linking)

The most typical use of the code in a C/C++ project is to use static linking. Having a sample code that uses SIM5

```C
#include <stdlib.h>
#include "sim5lib.h"

int main(int argc, char *argv[])
//***************************************************
{
    double a, rms, rbh;

    // for BH spin values from 0 to 1 print the radius of BH horizon and 
    // the location of the marginally stable orbit
    for (a=0.0; a<1.0; a+=0.01) {
        rbh = r_bh(a);  // radius of the event horizon
        rms = r_ms(a);  // radius of the marginally stable orbit
        fprintf(stderr, "%.4e  %e  %e\n", a, rbh, rms);
    }
    
    return 0;
}
```

we compile the code with (assuming SIM5 is located in subfolder sim5 of the current directory and assuming `make` has been run on SIM5)

```bash
gcc kerr-orbits.c -o kerr-orbits -I./sim5/lib -L./sim5/lib ./sim5/lib/sim5lib.o -lm
```

Refer to [examples](examples) for further inspiration.

### Using SIM5 in your Python project
SIM5 can be used as a Python module, which allows its routines to be called from a Python script. It is enough to `import sim5` in your script and call the member functions. Some care has to be taken when working with arrays and objects that are passed to SIM5 routines or returned from them. The SIM5 Python interface uses [SWIG](http://www.swig.org/) tool to make an interface between C and Python and you may refer to [SWIG documentation](http://www.swig.org/Doc3.0/SWIGDocumentation.html) to learn about how to do that.

A basic example of a Python script that calls SIM5 may look like this

```Python
import sys
import math
import sim5

bh_spin = 0.99    # BH spin parameter
bh_mass = 10.0    # BH mass in units of Solar mass
mdot    = 0.1     # disk accreation rate in units of Eddington accreation rate
alpha   = 0.1     # alpha parameter for viscosity
rout    = 500     # outer disk boudary in units of gravitational radius

# setup a relativistic thin disk model
sim5.disk_nt_setup(bh_mass, bh_spin, mdot, alpha, 0)

print '# Flux from the accretion disk surface as a function of radius for Novikov-Thorne disk model'
print '# mdot =', sim5.disk_nt_mdot(), '[Mdot_Edd]; ', sim5.disk_nt_mdot()*bh_mass*sim5.Mdot_Edd/1e18, 'g/s'
print '# format: radius [GM/c2] local flux [erg/s/cm2]'

r = sim5.disk_nt_r_min()
while (r < rout):
    print r, sim5.disk_nt_flux(r)
    r = r * 1.05
#end of while 
```
This example assumes that SIM5 subfolder does exist in the same folder along with the script. If that is not the case, put a line specifying the path to SIM5 library before the import statement.
```Python
sys.path.append('/path/to/sim5')
import sim5
...
```

Python calls to SIM5 routines are handled by the compiled C code in `lib/sim5lib.so` library. Although there is some overhead connected with the call to SIM5 routines, the actual evaluation of SIM5 functions is as fast as C code can be. This makes complex tasks (like raytracing) reasonably fast even if they are coded in Python (see examples in demo folder).

### Raytracing with SIM5

SIM5 library implements two essential ways of raytracing photons through curved spacetimes:
* by evaluating positional elliptic integrals along a defined trajectory
* by a stepwise integration of the [geodesic equation](https://en.wikipedia.org/wiki/Geodesics_in_general_relativity) for null geodesics

Both approaches have their advatages and disadvantages and the choice of one or the other depends mainly on the specific task that has to be carried out. In general, the use of the stepwise integration of the geodesic equation is necessary when the spacetime is not described by Kerr metric, otherwise the required accuracy, performance and computational complexity have to be considered to make a choice.

#### Raytracing by evaluating positional elliptic integrals

In Kerr spacetime, photon motion in *r*-*&theta;* plane is governed by equation
![equation of photon motion](https://latex.codecogs.com/gif.latex?\int^{r}\frac{dr}{\sqrt{R(r)}}=\int^{\theta}\frac{d\theta}{\sqrt{\Theta(\theta)}},)  
where *R(r)* and *&Theta;(&theta;)* are polynomials, and a solution can be expressed in terms of [Jacobi elliptic functions](https://en.wikipedia.org/wiki/Jacobi_elliptic_functions). Any photon geodesic is determined by it two constants of motion (traditionally denoted *&lambda;* and *q*) and a position along the geodesic can be uniquely expressed using the value of the above positional integral. In this way, once the trajectory is fixed by its motion constants, it is possible to follow the photon path by changing the value of the positional integral and by inverting it to obtain spacetime coordinates:
```C
// initialize the geodesic with an initial position `x` and direction `k`
geodesic_init_src(bh_spin, x[1], x[2], k, k[1]<0?1:0, &gd, &error);
// get the initial value of the positional integral
P = geodesic_P_int(&gd, x, 0);

r = x[1];   // radial coordinate
m = x[2];   // poloidal cordinate, m=cos(theta)
while (r < R_MAX) {
    // make a step along the geodesic; both `P`, `r` and `m` are updated
    geodesic_follow(&gd, deltaP, &P, &r, &m, &status);   
    if (!status) break;
    ...	 // do some stuff at the new place
}
```

Refer to the following examples:
* [tbd]


##### Equatorial plane
An important special case of the above approach is raytracing of photons that originate from the equatorial plane, as this is a commonly adopted approximation e.g. in case of dealing with accretion disks. Integrating over the observer's sky, we can directly determine the radius from which the photons originate by an analytic inversion of an elliptic positional integral. Thus, by using
```C
geodesic_init_inf(inclination, spin, alpha, beta, &gd, &error);
P = geodesic_find_midplane_crossing(&gd, 0);
r = geodesic_position_rad(&gd, P);
```
one can quickly and efficiently find the original radius of a photon coming out of the disk surface that has  impact parameters *alpha* and *beta* on the observer's image plane.

Refer to the following examples:
* [Example 04](examples/04-disk-image-eqplane) - accretion disk image
* [Example 05](examples/05-disk-spectrum-eqplane) - accretion disk spectrum

#### Raytracing by stepwise geodesic equation integration
Having an initial position *x<sup>&mu;</sup>* and direction of propagation *dx<sup>&mu;</sup>/d&lambda;* of the photon, its subsequent trajectory through the spacetime can be computed by solving a second-order ordinary differential equation (so called geodesic equation):
![geodesic equation](https://latex.codecogs.com/gif.latex?$%5Ccfrac%7Bd%5E2&space;x%5E%5Cmu%7D%7Bd%5Clambda%5E2%7D&space;=&space;-%5CGamma%5E%5Cmu_%7B%5Calpha%5Cbeta%7D%5Cbig%28x%5E%5Cmu%5Cbig%29%5C,%5Ccfrac%7Bdx%5E%5Calpha%7D%7Bd%5Clambda%7D&space;%5C,&space;%5Ccfrac%7Bdx%5E%5Cbeta%7D%7Bd%5Clambda%7D$)
SIM5 uses [Verlet integration](https://en.wikipedia.org/wiki/Verlet_integration) method to numericaly integrate geodesic equations because it offers good performance while keeping reasonable precission (see [Dolence at al. 2009](http://dx.doi.org/10.1088/0067-0049/184/2/387)). The basic scheme of the integration is the following (in Python):
```python
# initialize trajectory from initial position `x` and direction `k`
sim5.raytrace_prepare(bh_spin, x, k, 1.0, rtd)
# raytracing loop
while (True):
    dl.assign(1e-2)                    # set target step for affine parameter
    sim5.raytrace(x, k, dl, rtd)       # advance along photon trajectory
    if (rtd.error > 0.1): break        # check precission
    if (x[1] > MAX_DISTANCE): break    # check for terminating condition
#end of while
```

Refer to the following examples:
* [Example 04](examples/04) - accretion disk image - equatorial plane


### Code parallelization

As the library is a collection of routines, it does not have any parallelization built in itself, but it is very easy to write CPU/GPU parallel code with it. SIM5 uses strict encapsulation and all its functions use only the data that is passed in as parameters; there are no global variables in the library. This allows to write thread-safe codes that can be straightforwardly parallelized. Examples for OpenMP, MPI and CUDA parallelization are provided in the `examples` folder.

## Documentation

Documentation for the library is generated automatically (by [Doxygen](http://www.doxygen.nl/)) from the comments in the source code. See [doc/sim5lib-doc.md](doc/sim5lib-doc.md) for the compiled file with a list of routines and constants.

## Citing

If you want to refer to SIM5 in your paper, please give the following reference:
Bursa, M. 2017, Proceedings of Ragtime 17-19: Workshops on Black Holes and Neutron Stars, p.7 [[ADS link](http://adsabs.harvard.edu/abs/2017bhns.work....7B)]

If you have used SIM5 in your code or used a code that itself is using SIM5, please acknowledge the use of the software by [citing](http://ascl.net/wordpress/about-ascl/citing-the-ascl-and-codes/) its ASCL record: [ascl:1811.011](http://adsabs.harvard.edu/abs/2018ascl.soft11011B).



## License

SIM5 is released under the MIT License.

