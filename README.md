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

The compilation process creates a `lib` folder and puts two files there, `sim5lib.c` and `sim5lib.h`, which contain a complete C source code merged by the compiler process from all the individual source files in `src` folder. Use those two files to statically link SIM5 in your code (most typical usage in C/C++ codes).
<!--
The third file is a compiled library used for run-time linking. It is used by the Python interface, but you may use it for dynamical linking of SIM5 to your code too.
-->

### Compilation of the Python interface

The Python interface does have a few external dependecies and thus it does not compile by default when calling `make`. Before compiling the Python interface, make sure you have the following pre-requisities installed:

```bash
    # python heders
    apt install python-dev
    # SWIG - a tool that connects C/C++ programs with high-level programming languages
    apt install swig
```

Having those do

```bash
    make python
```

## Usage

The library can be compiled into a C/C++ code, and with a proper interface it can be used also in other languages(Fortran, ...). It also works as a Python package making most of its functions accessible from Python scripts with no extra effort. 

### Using SIM5 in your C/C++ project

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
        rbh = r_bh(a);
        rms = r_ms(a);
        fprintf(stderr, "%.4e  %e  %e\n", a, rbh, rms);
    }
    
    return 0;
}
```

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
This example assumes that SIM5 folder does exist in the same folder along with the script. If that is not the case, put a line specifying the path to SIM5 library before the import statement.
```Python
sys.path.append('/path/to/sim5')
import sim5
...
```

Python calls to SIM5 routines are handled by the compiled C code in `sim5lib.so` library. Although there is some overhead connected with the call to SIM5 routines, the actual evaluation of SIM5 functions is as fact as C code can be. This makes complex tasks (like raytracing) reasonably fast even if they are coded in Python (see examples in demo folder).

### Raytracing with SIM5

tbd

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

