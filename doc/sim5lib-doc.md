# SIM5 Library Reference

SIM5 is a collection of C routines for relativistic raytracing and radiation transfer. It has a special focus on raytracing from accretion disks, tori, hot spots or any custom 3D configuration of matter in Kerr geometry, but it can be used with any other metrics as well. It can handle both optically thick and thin sources as well as transport of polarization properties and helps to calculate the propagation of light rays from the source to an observer through a curved spacetimes. The library is threas-safe (with a few documented exceptions) compiles with Nvidia CUDA compiler which opens the door to massive parallelization using GPUs.

SIM5 also comes with a Python interface making it very easy to call its functions directly from Python scripts. In addition it provides few Python classes to handle some more complex tasks.

The following documentation provides a detailed reference to the functions of the library. The library also comes with couple of examples that illustrate how to piece the individual routines together to do some useful stuff. Routines are grouped together according to specific topics. The content gives a list of these topics with link to corresponding parts of the docs. If you are looking for something specific, use Ctrl-F.

<br>

## Content
* [Constants and unit conversions (sim5const)](#sim5const)
* [Thin disk routines (sim5disk-nt)](#sim5disk-nt)
* [Wrapper for dynamic linking an external library with a disk model (sim5disk)](#sim5disk)
* [Drawing numbers from statistical distributions (sim5distributions)](#sim5distributions)
* [Numerical integration (sim5integration)](#sim5integration)
* [Numerical interpolation (sim5interpolation)](#sim5interpolation)
* [Null geodesics in Kerr spacetime (sim5kerr-geod)](#sim5kerr-geod)
* [Basic spacetime routines (sim5kerr)](#sim5kerr)
* [Mathematical routines (sim5math)](#sim5math)
* [Mathematical macros (sim5math)](#sim5math)
* [Polarization properties of radiation (sim5polarization)](#sim5polarization)
* [Polynomial roots (sim5polyroots)](#sim5polyroots)
* [Radiative processes routines (sim5radiation)](#sim5radiation)
* [Raytracing (sim5raytrace)](#sim5raytrace)
* [Root finding (sim5roots)](#sim5roots)
* [Utility routines (sim5utils)](#sim5utils)


<br>

<br>

## <a name="sim5const"></a>sim5const.h - Constants and unit conversions



The library's default unit system is CGS, which is what astronomers like. Any unprefixed constant is given in this system of units. Sometimes constants in SI units (prefixed by `si_`) or in geometrical units (G=c=1; prefixed by `gu_`) are needed and few of them are defined too.

To convert between different units, some frequent conversion factors are defined (with the number "2" in the name). E.g. by multipying an energy value given in ergs by `erg2joule` factor, one gets the value in Joules. 


| Name | Description | Value |
|------|-------------|-------|
| grav_radius | gravitational radius GM/c2 of Sun [cm] | 1.476716e+05 |
| speed_of_light | speed of light [cm/s] | 2.997925e+10 |
| speed_of_light2 | square of speed of light [cm^2/s^2] | 8.987554e+20 |
| boltzmann_k | Boltzmann constant [erg/K] | 1.380650e-16 |
| sb_sigma | Stefan-Boltzmann constant [erg cm-2 s-1 K-4] | 5.670400e-05 |
| sigma_thomson | cross-section of Thomson scattering [cm^-2] | 6.652458e-25 |
| parsec | parsec [cm] | 3.085680e+18 |
| mass_proton | mass of proton [g] | 1.672622e-24 |
| mass_electron | mass of electron [g] | 9.109382e-28 |
| solar_mass | solar mass [g] | 1.988920e+33 |
| grav_const | gravitational constant [cm3 g-1 s-2] | 6.673000e-08 |
| planck_h | planck constant [erg.s] | 6.626069e-27 |
| atomic_mass_unit | atomic mass unit [g] | 1.660539e−24 |
| avogadro_number | Avogadro number [mol^-1] | 6.022141e+23 |
| Mdot_Edd | Eddington mass accretion rate [g/s * (M/Msun)] | 2.225475942e+18 |
| L_Edd | Eddington luminosity [erg/s * (M/Msun)] | 1.257142540e+38 |
| si_grav_radius | gravitational radius GM/c2 of Sun [m] | 1.476716e+03 |
| si_speed_of_light | speed of light [m/s] | 2.997925e+08 |
| si_speed_of_light2 | square of speed of light [m^2/s^2] | 8.987554e+16 |
| si_boltzmann_k | Boltzmann constant [J/K] | 1.380650e-23 |
| si_sb_sigma | Stefan-Boltzmann constant [J m-2 s-1 K-4] | 5.670400e-08 |
| si_electronvolt | electronvolt [J/elecronvolt] | 1.602177e-19 |
| si_parsec | parsec [m] | 3.085680e+16 |
| si_grav_const | gravitational constant [m3 kg-1 s-2] | 6.673000e-11 |
| si_erg | erg [J] | 1.000000e-07 |
| si_solar_mass | solar mass [kg] | 1.988920e+30 |
| si_angstrom | angstrom [meters] | 1.000000e-10 |
| si_sigma_T | Thomson scattering cross-section for electron [m^2] | 6.652459e-29 |
| si_mass_proton | mass of proton [kg] | 1.672622e-27 |
| si_mass_electron | mass of electron [kg] | 9.109382e-31 |
| si_planck_h | planck constant [J.s] | 6.626069e-34 |
| gu_sb_sigma | Stefan-Boltzmann constant (sb_sigma/(c^5/G)) [K-4] | 1.562539e-60 |
| erg2kev | erg->keV (erg/(1e3*electronvolt)) | 6.241507e+08 |
| kev2erg | keV-erg ((1e3*electronvolt)/erg) | 1.602177e-09 |
| joule2kev | Joule->keV (1/(1e3*electronvolt)) | 6.241507e+15 |
| joule2erg | Joule->erg | 1.000000e+07 |
| erg2joule | erg->Joule | 1.000000e-07 |
| kev2joule | keV->Joule (1e3*electronvolt) | 1.602177e-16 |
| freq2kev | Hz->keV (planck_h/(1e3*electronvolt)) | 4.135667e-18 |
| kev2freq | keV->Hz (1e3*electronvolt/planck_h) | 2.417990e+17 |
| msq2cmsq | squared centimeters in one squared meter | 1.000000e+04 |
| cmsq2msq | squared meters in one squared centimeter | 1.000000e-04 |
| kelvin2kev | Kelvin->keV (boltzmann_k/(1e3*electronvolt)) | 8.617342e-08 |
| kev2kelvin | keV->Kelvin ((1e3*electronvolt)/boltzmann_k) | 1.160451e+07 |
| m2cm | meter->centimeter | 1.000000e+02 |
| cm2m | centimeter->meter | 1.000000e-02 |
| kev2ev | keV->eV | 1.000000e+03 |
| ev2kev | eV->keV | 1.000000e-03 |
<br><br><br>

<br>

## <a name="sim5disk-nt"></a>sim5disk-nt.c - Thin disk routines



Provides routines for the radial structure of a relativistic thin disk model as given by Novikov & Thorne (1973) and Page & Thorne (1974).

NOTE: This unit uses some static variables to store some persistent information and due to that it is NOT thread-safe. For the same reasons, the routines declared here are not available to CUDA. 


#### disk_nt_setup()

Sets up a relativistic (Novikov-Thorne) model of a thin disk.


    DEVICEFUNC int disk_nt_setup(double M, double a, double mdot_or_L, double alpha, int options)

The disk can be set up using either mass accretion rate (mdot) or by specifying its luminosity (L). Mass accretion rate is passed as a ratio between the actual mass supply rate in grams per second to the Eddington mass acretion rate corresponding to the given mass M (the Eddington mass accretion rate and luminosity constant are declared in sim5const.h unit). In case of luminosity is the initial parameter, the model calculates the accretion rate in such a way that the integrated disk luminosity matches the given value in ergs/sec relative to the Eddington luminosity for the given black-hole mass.

**Parameters**

* **M**: mass of the central BH [M_sun] 
* **a**: spin of the central BH [0..1] 
* **mdot_or_L**: mass accretion rate (default) or luminosity (both in eddington units; see sim5const.h) 
* **alpha**: viscosity parameter 
* **options**: optional switches (default=0; switches can be combined with `+` operator)
  * DISK_NT_OPTION_LUMINOSITY: `mdot_or_L` parameter is interpreted as luminosity





**Return value**

A status code (currently returns always zero) 




<br>

---------

#### disk_nt_finish()

Finalize the disk model.


    DEVICEFUNC void disk_nt_finish()

Cleans up and frees necessary memory. 

<br>

---------

#### disk_nt_r_min()

Minimal radius of the disk (disk inner edge).


    DEVICEFUNC double disk_nt_r_min()

Provides minimal value for radius for which the functions provide valid results. For NT disk, this corresponds to the radius of the marginally stable orbit (r_ms, also known as ISCO), where there is zero torque in the fluid.

**Return value**

Radius of disk inner edge [GM/c2] 




<br>

---------

#### disk_nt_flux()

Local flux from one side of the disk.


    DEVICEFUNC double disk_nt_flux(double r)

Provides radial radiation flux dependence for Novikov-Thorne accretion disk. Formulae are based on Page&Thorne(1974) http://adsabs.harvard.edu/abs/1974ApJ...191..499P

Note the retuned flux is local flux, i.e. flux measured by an observer that is at rest with respect to the fluid.

**Parameters**

* **r**: radius of emission [GM/c2]


**Return value**

Total outgoing flux from unit area on one side of the disk [erg cm-2 s-1]. 




<br>

---------

#### disk_nt_lumi()

Total disk luminosity.


    DEVICEFUNC double disk_nt_lumi()

Luminosity is obtained by integrating local flux over the surface area of the disk (both sides) going into the whole sky (4pi solid angle). The integration makes a proper transformation of the flux from local to coordinate frame, but it ignores other relativistic effects, e.g. light bending.


```math
L = 2 * 2\pi \int F(r) (-U_t) r dr
```


**Return value**

Total disk luminosity of both surfaces [erg s-1] 




<br>

---------

#### disk_nt_mdot()

Mass accretion rate.


    DEVICEFUNC double disk_nt_mdot()

Returns mass accretion rate in Eddington units. See `disk_nt_setup()` for details.

**Return value**

Mass accretion rate in Eddington units. 




<br>

---------

#### disk_nt_sigma()

Column density.


    DEVICEFUNC double disk_nt_sigma(double r)

Returns midplane column density of the fluid, i.e. the fluid density integrated from midplane to the disk surface, at a given radius for the first two zones according to formulae from Chandrasekhar's book.


```math
 \Sigma = \int_0^H \rho dz 
```


**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Midplane column density in [g/cm2]. 




<br>

---------

#### disk_nt_ell()

Specific angular momentum.


    DEVICEFUNC double disk_nt_ell(double r)

Returns specific angular momentum of the fluid at given radius.

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Specific angular momentum in [g.u.]. 




<br>

---------

#### disk_nt_vr()

Radial velocity.


    DEVICEFUNC double disk_nt_vr(double r)

Returns bulk radial velocity of the fluid at given radius, which in case of thin disks is always zero.

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Radial velocity in [speed_of_light]. 




<br>

---------

#### disk_nt_h()

Surface height.


    DEVICEFUNC double disk_nt_h(double r)

Returns the scale-height of the surface of the disk (a measure of the effective photosphere location) above midplane at given radius. In thin disks, this is always zero. In fact, the height of the disk should be where the equation of hydrostatic equilibrium gives it, but the thin disk approximation assumes the disk razor thin, hence H=0.

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Scale-height [rg]. 




<br>

---------

#### disk_nt_dhdr()

Derivative of surface height.


    DEVICEFUNC double disk_nt_dhdr(double r)

Returns surface profile as derivative $`dH/dR`$ of its height above midplane at given radius. For thin disks, this is always zero.

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Derivative of surface height. 




<br>

---------

#### disk_nt_dump()

Prints the disk structure as a function of radius.


    DEVICEFUNC void disk_nt_dump(char *filename)

The function prints the profile of all quantities as a function of radius from r_ms to some outer radius (~2000 rg). It prints to a file identified by its path (overwrites existing) and if that is empty it prints to STDOUT.

**Parameters**

* **filename**: Path to a file that should be written. If NULL then it prints to STDOUT. 




<br>

---------

<br><br><br>

<br>

## <a name="sim5disk"></a>sim5disk.c - Wrapper for dynamic linking an external library with a disk model



Loads an external library and binds its functions. An external library (in Linux) is a compiled shared library with .so extension.

The external library has to declare at least the bellow listed methods and provide an implementation of them. Calls to those methods in SIM5 are passed to the linked library and the result is returned by a wrapper function. See a demo of how to use this functionality in the examples folder.

NOTE: This unit uses static variables to store persistent information about the linked library. As a result, routines in this module are NOT thread-safe in a sense different threads cannot each link a different library. They can, however, all make calls to the already linked library. For the same reasons, the routines declared here are not available to CUDA. 


#### diskmodel_init()

External disk model initialization.


    int diskmodel_init(char *modellib, double M, double a, char *params)

Loads the compiled library, links is into the program and calls its initialization routine. All the required functions have to be declared in the library.

**Parameters**

* **modellib**: filesystem path to the library (string) 
* **M**: mass of BH [Msun] 
* **a**: spin of BH [0..1] 
* **params**: parametres that are passed to the library initialization function


**Return value**

Return the result of the library's initialization function or -1 of the library could not loaded or be initialized. 




<br>

---------

#### diskmodel_done()

External disk model finitialization.


    void diskmodel_done()

Frees memory and unlinks the libraray. 

<br>

---------

#### diskmodel_name()

Model name.


    char * diskmodel_name()

Returns a pointer to a string with the model's name. 

<br>

---------

#### diskmodel_r_min()

Minimal radius of the disk (disk inner edge).


    double diskmodel_r_min()

Gives minimal value for radius for which the functions provide valid results. E.g. for NT disk, this corresponds to the radius of the marginally stable orbit.

**Return value**

Radius of disk inner edge [GM/c2] 




<br>

---------

#### diskmodel_mdot()

Mass accretion rate.


    double diskmodel_mdot()

Returns mass accretion rate in Eddington units of (Mdot_Edd*M).

**Return value**

Mass accretion rate in Eddington units. 




<br>

---------

#### diskmodel_lumi()

Total disk luminosity.


    double diskmodel_lumi()

Luminosity is obtained by integrating local flux over the surface area of the disk (both sides) going into the whole sky (4pi solid angle). The integration makes a proper transformation of the flux from local to coordinate frame, but it ignores other relativistic effects, e.g. light bending.


```math
L = 2 * 2\pi \int F(r) (-U_t) r dr
```


**Return value**

Total disk luminosity of both surfaces [erg s-1] 




<br>

---------

#### diskmodel_flux()

Local flux from one side of the disk.


    double diskmodel_flux(double R)

Provides radial radiation flux dependence measured in local frame, i.e. flux measured by an observer that is at rest with respect to the fluid.

**Parameters**

* **R**: radius of emission [GM/c2]


**Return value**

Total outgoing flux from unit area on one side of the disk [erg cm-2 s-1]. 




<br>

---------

#### diskmodel_sigma()

Column density.


    double diskmodel_sigma(double R)

Returns midplane column density of the fluid, i.e. the fluid density integrated from midplane to the disk surface, at a given radius.

**Parameters**

* **R**: radius (measured in equatorial plane) [rg]


**Return value**

Midplane column density in [g/cm2]. 




<br>

---------

#### diskmodel_ell()

Specific angular momentum.


    double diskmodel_ell(double R)

Returns specific angular momentum of the fluid at given radius.

**Parameters**

* **R**: radius (measured in equatorial plane) [rg]


**Return value**

Specific angular momentum in [g.u.]. 




<br>

---------

#### diskmodel_vr()

Radial velocity.


    double diskmodel_vr(double R)

Returns bulk radial velocity of the fluid at given radius as measured by aan observer in the co-rotating frame.

**Parameters**

* **R**: radius (measured in equatorial plane) [rg]


**Return value**

Radial velocity in [speed_of_light]. 




<br>

---------

#### diskmodel_h()

Surface height.


    double diskmodel_h(double R)

Returns the scale-height of the surface of the disk above midplane at given radius.

**Parameters**

* **R**: radius (measured in equatorial plane) [rg]


**Return value**

Scale-height [rg]. 




<br>

---------

#### diskmodel_dhdr()

Derivative of surface height.


    double diskmodel_dhdr(double R)

Returns surface profile as derivative $`dH/dR`$ of its height above midplane at given radius.

**Parameters**

* **R**: radius (measured in equatorial plane) [rg]


**Return value**

Derivative of surface height. 




<br>

---------

#### diskmodel_eval()

Other quantity evaluation.


    double diskmodel_eval(double R, int quantity)

Returns the value of a given quantity. The disk model may provide more quantities than the standard set. Additional quantities may be accessed using this function by providing the quantity identifier.

**Parameters**

* **R**: radius (measured in equatorial plane) [rg] 
* **quantity**: quantity identified code


**Return value**

The value of the requested quantity. 




<br>

---------

#### diskmodel_params()

Prints model parameters.


    void diskmodel_params(FILE *stream)/! if(stream) _diskmodel_params(stream)

Writes down the parameters of the model to the given stream. **Parameters**

* **stream**: stream to write to 




<br>

---------

#### diskmodel_dump()

Prints the disk structure as a function of radius.


    void diskmodel_dump(char *filename)

The function prints the profile of all quantities as a function of radius from r_ms to some outer radius (~2000 rg). It prints to a file identified by its path (overwrites existing) and if that is empty it prints to STDOUT.

**Parameters**

* **filename**: Path to a file that should be written. If NULL then it prints to STDOUT. 




<br>

---------

<br><br><br>

<br>

## <a name="sim5distributions"></a>sim5distributions.c - Drawing numbers from statistical distributions



Routines for drawing random numbers that follow a given statistical distribution. 


#### distrib_init()

Creates a distribution based on given probability density function (PDF).


    DEVICEFUNC void distrib_init(sim5distrib *d, double(*pdf)(double), double x_min, double x_max, int N)

Uses given probability density function to initiate the internal data structure with calculated cummulative distribution and its inverse.

**Parameters**

* **d**: pointer to structure that stores the disribution data 
* **pdf**: pointer to probability density function 
* **x_min**: left boundary for x 
* **x_min**: right boundary for x 
* **N**: number of samples to take over the interval [x_min:x_max] 




<br>

---------

#### distrib_done()

Frees the internal data for the distribution.


    DEVICEFUNC void distrib_done(sim5distrib *d)

**Parameters**

* **d**: pointer to structure that stores the disribution data 




<br>

---------

#### distrib_hit()

Generates a functional value on interval [x_min:x_max] according to the distribution.


    DEVICEFUNC INLINE double distrib_hit(sim5distrib *d)

With the help of precomputed cummulative distribution function it generates a value form the interval [x_min:x_max] and returns this value.

**Parameters**

* **d**: pointer to structure that stores the disribution data


**Return value**

value from the distribution 




<br>

---------

<br><br><br>

<br>

## <a name="sim5integration"></a>sim5integration.c - Numerical integration



Routines for simple numerical integrations of functions. 


#### integrate_trapezoid()

Integral of a function using trapezoid rule.


    DEVICEFUNC double integrate_trapezoid(double(*f)(double), double a, double b, double acc)

Computes the integral $` \int^a_b f(x) dx `$ in a series of refinement steps until the relative accuracy is better than requested or the maximum predefined number of steps is reached.

Note: Accuracy should not be increased beyond ~10^-6 as roundoff errors start to accumulate if too many steps are taken.

**Parameters**

* **f**: pointer fo a function that is to be integrated 
* **a**: interval of integration lower limit 
* **b**: interval of integration upper limit 
* **acc**: relative accuracy of integration


**Return value**

Integral of the function over the interval [a,b]. 




<br>

---------

#### integrate_simpson()

Integral of a function using Simpson rule.


    DEVICEFUNC double integrate_simpson(double(*f)(double), double a, double b, double acc)

Computes the integral $` \int^a_b f(x) dx `$ in a series of refinement steps until the relative accuracy is better than requested or the maximum predefined number of steps is reached. Simpson rule is generally more efficient than trapezoid rule when the function to be integrated has a finite 4th derivative (continuous 3rd derivative).

Note: Accuracy should not be increased beyond ~10^-6 as roundoff errors start to accumulate if too many steps are taken.

**Parameters**

* **f**: pointer fo a function that is to be integrated 
* **a**: interval of integration lower limit 
* **b**: interval of integration upper limit 
* **acc**: relative accuracy of integration


**Return value**

Integral of the function over the interval [a,b]. 




<br>

---------

#### gauleg()

Gauss-Legendre quadrature.


    DEVICEFUNC void gauleg(double x1, double x2, double x[], double w[], int n)

Given the lower and upper limits of integration $`x_1`$ and $`x_2`$, and given $`n`$ points, this routine returns arrays $`x[n]`$ and $`w[n]`$ of length $`n`$, containing the abscissas and weights of the Gauss-Legendre n-point quadrature formula.

To integrate a function $`f(x)`$ over the interval $`[x_1,x_2]`$, just do a simple summation: 
```
gauleg(x1, x2, x, w, n)
integral = sum(i=0, i<n-1, f(x[i])*w[i])
```



Using n points, the method exactly integrates polynomials up to (2*n-1)'th degree. If the function $`f(x)`$ is well approximated by a polynomial, the integral will be very accurate.

**Parameters**

* **x1**: interval of integration lower limit 
* **x2**: interval of integration upper limit 
* **x**: pointer to an array where abscissas should be written 
* **w**: pointer to an array where waights should be written 
* **n**: number of quandrature points and also the minimal dimensions of x[] and w[] arrays


**Return value**

Abscissas and weights are returned in x[] and w[] arrays. 




<br>

---------

#### qgaus()

Gauss-Legendre integral of a function.


    float qgaus(float(*func)(float), float a, float b)

Returns the integral of the function over the interval [a,b] computed using the by ten-point Gauss-Legendre integration: the function is evaluated exactly ten times at interior points in the range of integration.

**Parameters**

* **func**: pointer to a function to be integrated 
* **a**: interval of integration lower limit 
* **b**: interval of integration upper limit


**Return value**

Value of the integral over the interval [a,b]. 




<br>

---------

#### qgaus_general()

A generalized version of qgauss() function.


    double qgaus_general(double w[], double y[], int N, double a, double b)

It gives the integral of the function over the interval [a,b] on a supplied grid of functional values y and ther weights.

**Parameters**

* **func**: pointer to a function to be integrated 
* **w**: array of weights 
* **y**: array of functional values y=f(x[i]) 
* **N**: length of w[] and y[] arrays 
* **a**: interval of integration lower limit 
* **b**: interval of integration upper limit


**Return value**

Value of the integral over the interval [a,b]. 




<br>

---------

<br><br><br>

<br>

## <a name="sim5interpolation"></a>sim5interpolation.c - Numerical interpolation



Provides routines for interpolation of table data. 


#### sim5_interp_init()

Interpolation initialization.


    DEVICEFUNC void sim5_interp_init(sim5interp *interp, double xa[], double ya[], long N, int data_model, int interp_type, int interp_options)

Initializes the interpolation object `interp` with the data (xa,ya) where xa and ya are arrays of size N. **Parameters**

* **interp**: interpolation object 
* **xa**: array of x values (`xa` data array is always assumed to be strictly ordered, with increasing x values) 
* **ya**: array of y values 
* **N**: size of x/y arrays 
* **data_model**: switch that determines how input data in `xa` and `ya` arrays should be handled: INTERP_DATA_REF=X/Y arrays are referenced, INTERP_DATA_COPY=X/Y arrays are copied, INTERP_DATA_BUILD=X/Y arrays are not passed, they are build by calls to sim5_interp_data_push(); with INTERP_DATA_REF option (default), the interpolation object (interp) does not save the data arrays `xa` and `ya`, only saves pointers to them; with INTERP_DATA_COPY it makes an independent copy of those arrays, in which case the original arrays can be modified or freed after calling sim5_interp_init() 
* **interp_type**: determines in which way data will be interpolated: INTERP_TYPE_LINLIN=linear interpolation in both X and Y, INTERP_TYPE_LINLOG=linear interpolation in X, logarithmic in Y, INTERP_TYPE_LOGLIN=logarithmic interpolation in X, linear in Y, INTERP_TYPE_LOGLOG=logarithmic interpolation in both X and Y, INTERP_TYPE_SPLINE=linear cubic spline interpolation 
* **interp_options**: specifies additional options (a combination of options can be used): INTERP_OPT_ACCEL=interpolation will use acceleration (cashing of index values), INTERP_OPT_CAN_EXTRAPOLATE=extrapolation is allowed when an `x` value for an of-out-grid point is requested


**Return value**

Returns interp object to be used in actual interpolation. 




<br>

---------

#### sim5_interp_data_push()

Pushes data into interpolation object.


    DEVICEFUNC void sim5_interp_data_push(sim5interp *interp, double x, double y)

Adds a data point [x,y] into the interpolation object. This function is for filling interpolation object that has been created with option INTERP_DATA_BUILD with interolated data. The data must come in ordered sequence (x[i] < x[i+1])

**Parameters**

* **interp**: interpolation object 
* **x**: x-value of data point 
* **y**: y-value of data point 




<br>

---------

#### sim5_interp_eval()

Interpolated data evaluation.


    DEVICEFUNC double sim5_interp_eval(sim5interp *interp, double x)

Makes the evalutaion on interpolated grid at given point.

**Parameters**

* **interp**: interpolation object 
* **x**: value for which to get interpolated value


**Return value**

Interpolated value. 




<br>

---------

#### sim5_interp_done()

Interpolation finalization.


    DEVICEFUNC void sim5_interp_done(sim5interp *interp)

Frees the interpolation object interp (including a copied data, if necessary).

**Parameters**

* **interp**: interpolation object 




<br>

---------

#### sim5_interp_alloc()

Alloc interpolation object memory.


    DEVICEFUNC sim5interp * sim5_interp_alloc()

Makes a memory allocation for interpolation object. This function should be used for heap-allocated variant of usage: 
```
sim5interp* interp;
interp = sim5_interp_alloc();
sim5_interp_init(interp, ...);
sim5_interp_free(interp);
```



**Return value**

Interpolation object. 




<br>

---------

#### sim5_interp_free()

Free interpolation object memory.


    DEVICEFUNC void sim5_interp_free(sim5interp *interp)

Frees the interpolation object interp that had been previously alocated by `sim5_interp_alloc()`.

**Parameters**

* **interp**: interpolation object 




<br>

---------

<br><br><br>

<br>

## <a name="sim5kerr-geod"></a>sim5kerr-geod.c - Null geodesics in Kerr spacetime



Routines for computing null geodesics (trajectories of light rays) in Kerr spacetime based on numerical evaluation of elliptic integrals. 


#### geodesic_init_inf()

Initialization of a geodesic based on its impact parameters at infinity.


    DEVICEFUNC int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *error)

Makes a setup for a geodesic that is specified by its impact parameters at infinity. The imapct parameter is a perpendicular distance of the ray from the line-of-sight of the observer towards the black hole (the ray is parallel to the line-of-sight at infinity).

**Parameters**

* **i**: inclination angle of observer (angle between BH rotation axis and direction to observer) [radians] 
* **a**: BH spin [0..1] 
* **alpha**: impact parameter in horizontal direction [GM/c^2] 
* **beta**: impact parameter in vertical direction [GM/c^2] 
* **g**: structure with information about geodesic (output) 
* **error**: error code (output)


**Return value**

Returns TRUE on success or FALSE on error. In the latter case, an non-zero error code is returned in error parameter. Information about the geodesic is stored in structure g. 




<br>

---------

#### geodesic_init_src()

Initialization of a geodesic based on a position and direction.


    DEVICEFUNC int geodesic_init_src(double a, double r, double m, double k[4], int ppc, geodesic *g, int *error)

Makes a setup for a geodesic that is specified by a point and a direction (4-momentum vector) somewhere along the trajectory.

**Parameters**

* **a**: BH spin [0..1] 
* **r**: radial coordinate of the point [GM/c^2] 
* **m**: poloidal coordinate of the point ( $`m=cos(\theta)`$) 
* **k**: 4-momentum null vector ( $`k \cdot k=0`$) pointing in the direction of the ray 
* **ppc**: position with respect to pericenter (0=before pericenter, 1=after pericenter) 
* **g**: structure with information about geodesic (output) 
* **error**: error code (output)


**Return value**

Returns TRUE on success or FALSE on error. In the latter case, an non-zero error code is returned in error parameter. Information about the geodesic is stored in structure g. 




<br>

---------

#### geodesic_P_int()

Position integral along geodesics at radius r.


    DEVICEFUNC double geodesic_P_int(geodesic *g, double r, int ppc)

It gives the value of the integral (Bursa 2017, eq. 34, and 43) 
```math
 P = \int 1/\sqrt{R} dr = \int 1/\sqrt{\Theta} d\theta 
```
 The integral is integrated from infinity to the given point on the trajecotry, where the geodesic reaches radius $`r`$ either before or behind the trajecory pericenter.

The value of the integral increases monotonicly from infinity along the geodesic and thus it provides a convenient way of parametrizing the position along the geodesic that is used in many other routined of the module. Note however, that the value of this integral is not the affine parameter, which would be another choice for parametrization.

**Parameters**

* **g**: geodesic data 
* **r**: radial coordinate [GM/c^2] 
* **ppc**: position with respect to pericenter (0=before pericenter, 1=after pericenter)


**Return value**

Value of the position integral between infinity and given point. 




<br>

---------

#### geodesic_position_rad()

Radius at which the position integral gains value P.


    DEVICEFUNC double geodesic_position_rad(geodesic *g, double P)

Given the value $`P`$ of the positional integral along the geodesic, the function computes the radial coordinate for the position.

**Parameters**

* **g**: geodesic data 
* **P**: value of the position integral


**Return value**

Radial coordinate value [GM/c^2] or NAN in case of error. 




<br>

---------

#### geodesic_position_pol()

Poloidal coordinate value at which the position integral gains value P.


    DEVICEFUNC double geodesic_position_pol(geodesic *g, double P)

Given the value $`P`$ of the positional integral along the geodesic, the function computes the poloidal coordinate for the position. The coordinate is returned as a cosine of the angle theta.

**Parameters**

* **g**: geodesic data 
* **P**: value of the position integral


**Return value**

Poloidal coordinate value [cos(theta)] or NAN in case of error. 




<br>

---------

#### geodesic_position_pol_sign_k_theta()

Sign of the $k^\theta$ component of the 4-momentum.


    DEVICEFUNC double geodesic_position_pol_sign_k_theta(geodesic *g, double P)

Gives the orientation of the 4-momentum vector in the poloidal direction by returning the sign of $`k^\theta`$ component of the momentum vector.

**Parameters**

* **g**: geodesic data 
* **P**: value of the position integral


**Return value**

Returns +1 or -1 or NAN in case of an error. 




<br>

---------

#### geodesic_position_azm()

Azimuthal coordinate value at which the position integral gains value P.


    DEVICEFUNC double geodesic_position_azm(geodesic *g, double r, double m, double P)

Given the value $`P`$ of the positional integral along the geodesic, the function computes the azimuthal coordinate for the position. The value of azimuthal angle is assumed to be zero at infinity and the function gives the change of the angle between the point [r,m] and infinity.

**Parameters**

* **g**: geodesic data 
* **r**: radial coordinate of the current position 
* **m**: poloidal coordinate of the current position 
* **P**: value of the position integral


**Return value**

Difference between azimuthal angle at [r,m] and at infinity [radians]. The result can be larger than 2pi. 




<br>

---------

#### geodesic_timedelay()

Time delay (travel time) between positions P1 and P2.


    DEVICEFUNC double geodesic_timedelay(geodesic *g, double P1, double r1, double m1, double P2, double r2, double m2)

Gives time it takes the light to travel between two points along a geodesic. Returned value is always positive independent of the relative position of P1 and P2 along the geodesic.

**Parameters**

* **g**: geodesic data 
* **P1**: value of the position integral at point A 
* **r1**: value of the radial coordinate at point A; if zero, it is computed from P1 internally 
* **m1**: value of the poloidal coordinate at point A ( $`m=cos(\theta)`$); if zero, it is computed from P1 internally 
* **P2**: value of the position integral at point B 
* **r2**: value of the radial coordinate at point B; if zero, it is computed from P2 internally 
* **m2**: value of the poloidal coordinate at point B ( $`m=cos(\theta)`$); if zero, it is computed from P2 internally


**Return value**

Timedelay (positive) between position P1 and P2 (point A and B) or NAN in case of an error. 




<br>

---------

#### geodesic_dm_sign()

Gives the sign of the derivative d(m)/d(P) at current position.


    DEVICEFUNC double geodesic_dm_sign(geodesic *g, double P)

Parameters: g - geodesic P - value of the position integral

Returns: Sign of d(m)/d(P), i.e. +1 or -1. 

<br>

---------

#### geodesic_momentum()

Photon 4-momentum.


    DEVICEFUNC void geodesic_momentum(geodesic *g, double P, double r, double m, double k[])

Gives the 4-momentum of photons at given position along the geodesic. The function needs to know [r,m] coordinates of the point at the trajectory. If both r=m=0.0, the required values are computed from the value P of the position integral. To save these computations, value of [r,m] coordinates can be given to the function, if they have been computed before. Note: It is important to give the correct values of [r,m] corresponding to current position. The orientation of the momentum vector is always in the direction of increasing P, i.e. it points towards the radial turning point before it is reached and away from the radial turning point after it is reached.

**Parameters**

* **g**: - geodesic data 
* **P**: - value of the position integral 
* **r**: - radial coordinate (value or zero) 
* **m**: - cosine of poloidal coordinate (value or zero) 
* **k**: - 4-momentum vector (output)


**Return value**

Photon 4-momentum vector is returned in k[]. 




<br>

---------

#### geodesic_find_midplane_crossing()

Finds a crossing of the geodesic with the equatorial plane.


    DEVICEFUNC double geodesic_find_midplane_crossing(geodesic *g, int order)

Calculates, where (if ever) the geodesic crosses the equatorial plane, and returns the value of the positional integral for that place. This is the fastest way to integrate over the equatorial plane. The positional integral can be converted to radius, which allows straightforward integration over the solid angle, for example 
```math
 F_\nu(E) = 1/D^2 \int I_\nu(E/g, r) d\alpha\,d\beta 
```
 where $`r = r(\alpha, \beta)`$, $`D`$ is distance and $`g`$ is g-factor.

**Parameters**

* **g**: geodesic data 
* **order**: order of crossing; order=0 zero is the first crossing, higher orders may be reached by some geodesics that loop around the photon orbit


**Return value**

Value of the positional integral at the equatorial plane. 




<br>

---------

#### geodesic_follow()

Makes a step along the geodesic.


    DEVICEFUNC void geodesic_follow(geodesic *g, double step, double *P, double *r, double *m, int *status)

Moves current position on the ray along the geodescis. On input, the function receives the current position integral value and radial and poloidal coordinate. It computes new values for P, r and m and shifts them to a new position along the geodesics by a given step.

This function is meant to be called in a cycle to follow the geodesics in a piecewise steps.

**Parameters**

* **g**: geodesic data 
* **step**: size of step to advance 
* **P**: value of the positional integral (input and output) 
* **r**: value of the radial coordinate (input and output) 
* **m**: value of the poloidal coordinate (input and output; $`m=cos(theta)`$) 
* **status**: status code; status=0 if ok, it get a non-zero value on an error 




<br>

---------

<br><br><br>

<br>

## <a name="sim5kerr"></a>sim5kerr.c - Basic spacetime routines



Provides basic routines for doing physics in Kerr and Minkowski spacetimes (metric, connection, tetrads, vector algebra, orbital motion, photon propagation, etc). 


#### flat_metric()

Flat (Minkowski) spacetime metric.


    DEVICEFUNC void flat_metric(double r, double m, sim5metric *metric)

Returns covariant Minkowski metric $`\eta_\mu\nu`$ in spherical coordinates (t,r,theta,phi).

**Parameters**

* **r**: radial coordinate 
* **m**: poloidal coordinate $`m=\cos\theta`$


**Return value**

Metric components are returned in metric parameter. 




<br>

---------

#### flat_metric_contravariant()

Flat (Minkowski) spacetime metric (contravariant).


    DEVICEFUNC void flat_metric_contravariant(double r, double m, sim5metric *metric)

Returns contravariant Minkowski metric $`\eta^\mu\nu`$ in spherical coordinates (t,r,theta,phi).

**Parameters**

* **r**: radial coordinate 
* **m**: poloidal coordinate $`m=\cos\theta`$


**Return value**

Metric components are returned in metric parameter. 




<br>

---------

#### kerr_metric()

Kerr spacetime metric.


    DEVICEFUNC void kerr_metric(double a, double r, double m, sim5metric *metric)

Returns Kerr metric $`g_\mu\nu`$ in spherical coordinates (t,r,theta,phi).

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate 
* **m**: poloidal coordinate $`m=\cos\theta`$


**Return value**

Metric components are returned in metric parameter. 




<br>

---------

#### kerr_metric_contravariant()

Kerr spacetime metric (contravariant).


    DEVICEFUNC void kerr_metric_contravariant(double a, double r, double m, sim5metric *metric)

Returns contravariant Kerr metric $`g^\mu\nu`$ in spherical coordinates (t,r,theta,phi).

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate 
* **m**: poloidal coordinate $`m=\cos\theta`$


**Return value**

Metric components are returned in metric parameter. 




<br>

---------

#### kerr_connection()

Christoffel symbol components for Kerr metric ( $Gamma^\mu_\alpha\beta$).


    DEVICEFUNC void kerr_connection(double a, double r, double m, double G[4][4][4])

Returns a matrix of connection coefficients for Kerr metric.

(!) NOTE: Christoffel tensor Gamma^i_(jk) is symmetric in lower two indices (j,k). This function only evaluates components of the tensor, where j<k and multiplies the value of these coeficients by a factor of two. In this way, both summation over j=1..4,k=1..4 and over j=1..4,k=j..4 will give the same result; however, care must be taken when summing Gamma^i_(jk) U^j V^k where factor 0.5 has to be used: `0.5*G[i][j][k]*(u[j]*v[k] + u[k]*v[j])`

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate 
* **m**: poloidal coordinate $`m=\cos\theta`$


**Return value**

Connection coeficients are returned in G parameter. 




<br>

---------

#### flat_connection()

Christoffel symbol components for Minkowski metric ( $Gamma^\mu_\alpha\beta$).


    DEVICEFUNC void flat_connection(double r, double m, double G[4][4][4])

Returns a matrix of connection coefficients for Minkowski metric, i.e. Kerr metric in the limit of M=0 and a=0.

(!) NOTE: Christoffel tensor Gamma^i_(jk) is symmetric in lower two indices (j,k). This function only evaluates components of the tensor, where j<k and multiplies the value of these coeficients by a factor of two. In this way, both summation over j=1..4,k=1..4 and over j=1..4,k=j..4 will give the same result; however, care must be taken when summing Gamma^i_(jk) U^j V^k where factor 0.5 has to be used: `0.5*G[i][j][k]*(u[j]*v[k] + u[k]*v[j])`

**Parameters**

* **r**: radial coordinate 
* **m**: poloidal coordinate $`m=\cos\theta`$


**Return value**

Connection coeficients are returned in G parameter. 




<br>

---------

#### Gamma()

Change of vector along a trajectory.


    DEVICEFUNC INLINE void Gamma(double G[4][4][4], double U[4], double V[4], double result[4])

Returns product of summation `-G^i_(jk) U^j V^k`. This is useful for calculating parallel transport of vector `V` along a trajectory specified by tangent vector `U`, as it gives the derivative $`dV/d\lambda`$. (Note the definition of G tensor components in symmetric indices.)

**Parameters**

* **G**: connection coefficients 
* **U**: tangent vector 
* **V**: transported vector


**Return value**

Change of transported vector $dV/d\lambda$. 




<br>

---------

#### vector_set()

Make a 4-vector.


    DEVICEFUNC INLINE void vector_set(double x[4], double x0, double x1, double x2, double x3)

Returns a 4-vector that has given components (contravarient form $`X^\mu`$).

**Parameters**

* **x0-x3**: components of the vector


**Return value**

Vector is returned in x parameter. 




<br>

---------

#### vector_copy()

Copy a 4-vector.


    DEVICEFUNC INLINE void vector_copy(double src[4], double dst[4])

Copies a 4-vector `src` into `dst`.

**Parameters**

* **src**: the source vector 
* **dst**: the target vector


**Return value**

A copy of src vector in dst. 




<br>

---------

#### vector_covariant()

Covariant version of a vector.


    DEVICEFUNC INLINE void vector_covariant(double V1[4], double V2[4], sim5metric *m)

Converts a standard (contravariant) form of a vector to covariant form ( $`X^\mu -> X_\mu`$). Metric `m` can be NULL in which case flat metric is assumed.

**Parameters**

* **V1**: contravariant vector 
* **V2**: covariant vector (output) 
* **m**: metric


**Return value**

Transformed vector is returned in V2 parameter. 




<br>

---------

#### vector_norm()

Norm of a vector.


    DEVICEFUNC INLINE double vector_norm(double V[4], sim5metric *m)

Returns norm of a vector V, i.e. sqrt(V*V). Metric `m` can be NULL in which case flat metric is assumed. NOTE: only works for space-like vectors, where V*V>0.

**Parameters**

* **V**: vector 
* **m**: metric


**Return value**

Norm of a vector. 




<br>

---------

#### vector_3norm()

3-norm of a vector (in flat spacetime).


    DEVICEFUNC INLINE double vector_3norm(double V[4])

Returns a 3-norm of a vector V (includes only 'spatial' components) in flat (Minkowski) spacetime, i.e. sqrt( V[i]*V[i]), where i=1..3.

**Parameters**

* **V**: vector 
* **m**: metric


**Return value**

3-norm of a vector. 




<br>

---------

#### vector_multiply()

Multiply a vector.


    DEVICEFUNC INLINE void vector_multiply(double V[4], double factor)

Multiplies the components of a vector by given factor.

**Parameters**

* **V**: vector 
* **factor**: multiplication factor


**Return value**

Modified vector V. 




<br>

---------

#### vector_norm_to()

Normalizes a vector to given size.


    DEVICEFUNC INLINE void vector_norm_to(double V[4], double norm, sim5metric *m)

Changes the vector components such that the norm of the vector becomes `|norm|`, i.e. $`V.V = \rm{norm}`$. The vector `V` must be space-like if norm>0, and it must be time-like if norm<0 (no checks are performed to enfoce that). Metric `m` can be NULL in which case flat metric is assumed. NOTE: vector cannot be null vector (V*V=0), nor norm can be zero; use vector_norm_to_null() in that case.

**Parameters**

* **V**: vector to be normalized 
* **norm**: normalization the vector should have 
* **m**: metric


**Return value**

Changes the components of vector V. 




<br>

---------

#### vector_norm_to_null()

Normalizes a null vector to given size.


    DEVICEFUNC void vector_norm_to_null(double V[4], double V0, sim5metric *m)

Changes the vector components such that the time-component is set to V0 and the remaining components are adjusted such the vector is still a null vector (V.V=0). Metric `m` can be NULL in which case flat metric is assumed. NOTE: the original vector must be null vector (V*V=0) and `V0` cannot be zero (no check is performed on that).

**Parameters**

* **V**: vector to be normalized 
* **V0**: the new value of vector time-component 
* **m**: metric


**Return value**

Changes the components of vector V. 




<br>

---------

#### dotprod()

Scalar product.


    DEVICEFUNC INLINE double dotprod(double V1[4], double V2[4], sim5metric *m)

Calculates scalar product of two 4-vectors ( $`U^\mu V^\nu g_{\mu\nu}`$). NOTE: if metric is NULL, flat space metric is used

**Parameters**

* **V1**: first vector 
* **V1**: second vector 
* **m**: metric


**Return value**

The value of scalar product 




<br>

---------

#### tetrad_general()

Tetrad of an observer that is moving with a general 4-velocity U.


    DEVICEFUNC void tetrad_general(sim5metric *m, double U[], sim5tetrad *t)

Note: the orientation of theta-vector points in opposite direction than d/d(theta), i.e. at the equatorial plane the theta-vector points in upward direction. (based on Kulkarni+2011, Dexter2016 eq.36-43)

**Parameters**

* **m**: metric 
* **U**: observer 4-velocity 
* **t**: tetrad


**Return value**

Returns tetrad vectors in t. 




<br>

---------

#### tetrad_zamo()

Tetrad of a zero angular momentum observer (ZAMO).


    DEVICEFUNC void tetrad_zamo(sim5metric *m, sim5tetrad *t)

Returns basis vectors for ZAMO (locally non-rotating) observer tetrad e_{(a)}^{}. Note the orientation of theta-vector that points in opposite direction than d/d(theta), i.e. at the equatorial plane the theta-vector points in upward direction.

**Parameters**

* **m**: metric 
* **t**: tetrad


**Return value**

Returns tetrad vectors in t. 




<br>

---------

#### tetrad_radial()

Tetrad of observer that moves purely in radial direction.


    DEVICEFUNC void tetrad_radial(sim5metric *m, double v_r, sim5tetrad *t)

Returns basis vectors of a tetrad e_{(a)}^{} of an observer that moves purely in radial direction (it does have zero phi velocity component). Positive value of `v_r` means outward motion, a negative value inward motion.

Tetrad vector orientation: e[1] (x-vector) is oriented along increasing r (even when off equatorial plane), e[2] (z-vector) is oriented along decreasing theta (goes "upwards" from eq plane), e[3] (y-vector) is oriented along increasing phi.

**Parameters**

* **m**: metric 
* **v_r**: radial velocity component in [c] 
* **t**: tetrad


**Return value**

Returns tetrad vectors in t. 




<br>

---------

#### tetrad_azimuthal()

Tetrad of observer that moves purely in azimuthal direction.


    DEVICEFUNC void tetrad_azimuthal(sim5metric *m, double Omega, sim5tetrad *t)

Returns basis vectors of a tetrad e_{(a)}^{} of an observer that moves purely in azimuthal direction with angular velocity Omega.

Tetrad vector orientation: e[1] (x-vector) is oriented along increasing r (even when off equatorial plane), e[2] (z-vector) is oriented along decreasing theta (goes "upwards" from eq plane), e[3] (y-vector) is oriented along increasing phi.

**Parameters**

* **m**: metric 
* **Omega**: angular velocity in [g.u.] 
* **t**: tetrad


**Return value**

Returns tetrad vectors in t. 




<br>

---------

#### tetrad_surface()

Tetrad of observer that moves along a surface.


    DEVICEFUNC void tetrad_surface(sim5metric *m, double Omega, double V, double dhdr, sim5tetrad *t)

Returns basis vectors of a tetrad e_{(a)}^{} of an observer that moves along an axisymmetric surface in Kerr spacetime, i.e. the obeserver moves azimuthally with angular velocity Omega and drifs radially along the surface with velocity V, which itself is measured locally in the corotating frame. The orientation of the surface is given locally by the derivative dH/dR, where H is the height of the surface above equatoriual plane [r*cos(theta)] and R is the BL radial coordinate in equatorial plane [r*sin(theta)].

Based on: Sadowski+2011, Appendix A (http://adsabs.harvard.edu/abs/2011A%26A...532A..41S); note the + metric signature used there.

Tetrad vector orientation: e[0] (t-vector) is oriented along the direction of the observer's 4-velocity e[1] (x-vector) is a spacelike vector in [r,theta] plane tangent to the surface and oriented outwards e[2] (z-vector) is a spacelike vector in [r,theta] plane normal to the surface and oriented upwards e[3] (y-vector) is a remaining spacelike vector orthogonal to all the others

**Parameters**

* **m**: metric 
* **Omega**: angular velocity in [g.u.] 
* **radial**: drift velocity mesured in corotating frame [c] 
* **t**: tetrad


**Return value**

Returns tetrad vectors in t. 




<br>

---------

#### bl2on()

Vector transformation from coordinate to local frame.


    DEVICEFUNC void bl2on(double Vin[4], double Vout[4], sim5tetrad *t)

Transforms a vector `V` from coordinate (Boyer-Lindquist) frame to local (orthonormal) frame specified by tetrad `t`.

MATH: V^(a) = e^(a)_ * V^, where e^(a)_ = e^(b) * g_ * n^ab V^(a) = dotprod(e_(b)^, Vin^) * n^ab

**Parameters**

* **Vin**: vector to transform (in coordinate basis) 
* **Vout**: transformed vector (in local basis) 
* **t**: transformation tetrad


**Return value**

Local vector Vout. 




<br>

---------

#### on2bl()

Vector transformation from local to coordinate frame.


    DEVICEFUNC void on2bl(double Vin[4], double Vout[4], sim5tetrad *t)

Transforms a vector `V` from local (orthonormal) frame specified by tetrad `t` to coordinate (Boyer-Lindquist) frame.

MATH: V^ = e^(a) * V^(a) or V^i = V^j * e_j^i

**Parameters**

* **Vin**: vector to transform (in local basis) 
* **Vout**: transformed vector (in coordinate basis) 
* **t**: transformation tetrad


**Return value**

Coordinate vector Vout. 




<br>

---------

#### r_bh()

Black hole event horizon radius.


    DEVICEFUNC INLINE double r_bh(double a)

**Parameters**

* **a**: black hole spin


**Return value**

Event horizon radius in [rg]. 




<br>

---------

#### r_ms()

Radius of marginally stable orbit.


    DEVICEFUNC INLINE double r_ms(double a)

**Parameters**

* **a**: black hole spin


**Return value**

Marginally stable orbit radius in [rg]. 




<br>

---------

#### r_mb()

Radius of marginally bound orbit.


    DEVICEFUNC INLINE double r_mb(double a)

Marginally bound orbit is minimal radius of bound (E<0) circular and parabolic orbits. http://adsabs.harvard.edu/abs/1972ApJ...178..347B, eq. 2.19.

**Parameters**

* **a**: black hole spin


**Return value**

Marginally bound orbit radius in [rg]. 




<br>

---------

#### r_ph()

Radius of photon orbit.


    DEVICEFUNC INLINE double r_ph(double a)

Photon orbit is radius of unstable circular photon orbit. http://adsabs.harvard.edu/abs/1972ApJ...178..347B, eq. 2.18.

**Parameters**

* **a**: black hole spin


**Return value**

Marginally bound orbit radius in [rg]. 




<br>

---------

#### OmegaK()

Angular frequency of Keplerian orbital motion.


    DEVICEFUNC INLINE double OmegaK(double r, double a)

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin


**Return value**

Angular frequency [Hz]. 




<br>

---------

#### ellK()

Specific angular momentum of Keplerian orbital motion.


    DEVICEFUNC INLINE double ellK(double r, double a)

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin


**Return value**

Keplerian specific angular momentum [g.u.]. 




<br>

---------

#### omega_r()

Angular frequency of radial epyciclic motion.


    DEVICEFUNC INLINE double omega_r(double r, double a)

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin


**Return value**

Angular velocity [Hz]. 




<br>

---------

#### omega_z()

Angular frequency of vertical epyciclic motion.


    DEVICEFUNC INLINE double omega_z(double r, double a)

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin


**Return value**

Angular velocity [Hz]. 




<br>

---------

#### Omega_from_ell()

Angular frequency corresponding to given angular momentum.


    DEVICEFUNC INLINE double Omega_from_ell(double ell, sim5metric *m)

**Parameters**

* **ell**: specific angular momentum [g.u.] 
* **m**: metric


**Return value**

Angular frequency [Hz]. 




<br>

---------

#### ell_from_Omega()

Specific angular momentum corresponding to given angular frequency.


    DEVICEFUNC INLINE double ell_from_Omega(double Omega, sim5metric *m)

**Parameters**

* **Omega**: angular frequency [Hz] 
* **m**: metric


**Return value**

Specific angular momentum [g.u.]. 




<br>

---------

#### gfactorK()

Redshift factor.


    DEVICEFUNC INLINE double gfactorK(double r, double a, double l)

Relativistic correction to energy of photon that is emitted by fluid at Keplerian rotation in equatorial plane; includes Doppler effect and gravitational redshift.

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin 
* **l**: photon motion constant lambda


**Return value**

Redshift factor. 




<br>

---------

#### photon_momentum()

Photon 4-momentum vector .


    DEVICEFUNC void photon_momentum(double a, double r, double m, double l, double q, double r_sign, double m_sign, double k[4])

Returns photon 4-momentum vector k^ such that k*k=0 (null vector). The orientation of the resulting vector `k` is given by the signs of `r_sign` and `m_sign` parameters. See MTW; Rauch & Blandford (1994), Appendix A; Li+05

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate [rg] 
* **m**: poloidal coordinate [cos(theta)] 
* **l**: photon motion constant lambda 
* **q**: photon motion constant Q (Carter's constant) 
* **r_sign**: sign of k[1] component of resulting momentum vector 
* **m_sign**: sign of k[2] component of resulting momentum vector 
* **k**: resulting momentum vector (output)


**Return value**

Photon momentum vector k. 




<br>

---------

#### photon_motion_constants()

Constants of motion L,Q for null geodesic.


    DEVICEFUNC void photon_motion_constants(double a, double r, double m, double k[4], double *L, double *Q)

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate [rg] 
* **m**: poloidal coordinate [cos(theta)] 
* **k**: photon 4-momentum vector 
* **L**: photon motion constant lambda (output) 
* **Q**: photon motion constant Q^2 (Carter's constant, output)


**Return value**

Photon motion constants L and Q. 




<br>

---------

#### photon_carter_const()

Carter's constant Q for null geodesic.


    DEVICEFUNC double photon_carter_const(double k[4], sim5metric *metric)

**Parameters**

* **k**: photon 4-momentum vector 
* **m**: metric


**Return value**

Carter's constants Q. 




<br>

---------

#### fourvelocity_zamo()

4-velocity of ZAMO (locally non-rotating) observer.


    DEVICEFUNC INLINE void fourvelocity_zamo(sim5metric *m, double U[4])

**Parameters**

* **m**: metric 
* **U**: 4-velocity (output)


**Return value**

4-velocity U. 




<br>

---------

#### fourvelocity_azimuthal()

4-velocity of azimuthally rotating observer.


    DEVICEFUNC INLINE void fourvelocity_azimuthal(double Omega, sim5metric *m, double U[4])

**Parameters**

* **Omega**: angular frequency of circular motion 
* **m**: metric 
* **U**: 4-velocity (output)


**Return value**

4-velocity U. 




<br>

---------

#### fourvelocity_radial()

4-velocity of radially moving observer.


    DEVICEFUNC INLINE void fourvelocity_radial(double vr, sim5metric *m, double U[4])

**Parameters**

* **vr**: radial component of 4-velocity (dr/dtau) 
* **m**: metric 
* **U**: 4-velocity (output)


**Return value**

4-velocity U. 




<br>

---------

#### fourvelocity_norm()


    DEVICEFUNC INLINE double fourvelocity_norm(double U1, double U2, double U3, sim5metric *m)

<br>

---------

#### fourvelocity()


    DEVICEFUNC void fourvelocity(double U1, double U2, double U3, sim5metric *m, double U[])

<br>

---------

<br><br><br>

<br>

## <a name="sim5math"></a>sim5math.c - Mathematical routines




#### sim5round()


    DEVICEFUNC INLINE long sim5round(double num)

<br>

---------

#### factorial()


    DEVICEFUNC INLINE long int factorial(long int n)

<br>

---------

#### reduce_angle_pi()


    DEVICEFUNC INLINE double reduce_angle_pi(double phi)

<br>

---------

#### reduce_angle_2pi()


    DEVICEFUNC INLINE double reduce_angle_2pi(double phi)

<br>

---------

#### ensure_range()


    DEVICEFUNC int ensure_range(double *val, double min, double max, double acc)

<br>

---------

#### sim5seed()


    DEVICEFUNC INLINE void sim5seed()

<br>

---------

#### sim5rand()


    DEVICEFUNC INLINE unsigned long long sim5rand()

<br>

---------

#### sim5urand()


    DEVICEFUNC INLINE double sim5urand()

<br>

---------

#### cartesian2spherical1()


    DEVICEFUNC void cartesian2spherical1(double x, double y, double z, double Vx, double Vy, double Vz, double *Vr, double *Vh, double *Vf)

<br>

---------

#### cartesian2spherical2()


    DEVICEFUNC void cartesian2spherical2(double cos_h, double sin_f, double cos_f, double Vx, double Vy, double Vz, double *Vr, double *Vh, double *Vf)

<br>

---------

#### makeComplex()


    DEVICEFUNC INLINE sim5complex makeComplex(double r, double i)

<br>

---------

#### nullComplex()


    DEVICEFUNC INLINE sim5complex nullComplex()

<br>

---------

<br><br><br>

<br>

## <a name="sim5math"></a>sim5math.h - Mathematical macros



Some useful mathematical macros. 


| Name | Description | Value |
|------|-------------|-------|
| PI | PI | 3.14159265359 |
| PI2 | 2*PI | 6.28318530718 |
| PI4 | 4*PI | 12.5663706144 |
| PI_half | PI/2 | 1.57079632679 |
| sqr | quadratic power of a | ((a) * (a)) |
| sqr2 | quadratic power of a | ((a) * (a)) |
| sqr3 | cubic power of a | ((a) * (a) * (a)) |
| sqr4 | quartic power of a | ((a) * (a) * (a) * (a)) |
| sqrt3 | cubic root of a | cbrt(a) |
| sqrt4 | quartic root of a | pow(a,0.25) |
| max | maximum of two values | ((a) > (b) ? (a) : (b)) |
| min | minimum of two values | ((a) < (b) ? (a) : (b)) |
| minmax | values within limits | min(vmax,max(val,vmin)) |
| odd | is odd number | ((a%2==1)?1:0) |
| sign | positive or negative | ((a) >= 0.0 ? (+1.0) : (-1.0)) |
| deg2rad | convert degrees to radians | ((a)/180.0*M_PI) |
| rad2deg | convert radians to degrees | ((a)*180.0/M_PI) |
| EE | power of 10 | pow(10.0,a) |
| ave | weighted average of two values | ((1.0-(w))*(a) + (w)*(b)) |
| logave | weighted logarithmic average of two values | (exp((1.0-(w))*log(a) + (w)*log(b))) |
| inrange | value in range | (((a)>=(min))&&((a)<=(max))) |
| rand | integer random number (long long int) | sim5rand() |
| urand | random number from 0 to 1 (double) | sim5urand() |
| ComplexI | complex unit (sqrt(-1)) | _Complex_I |
<br><br><br>

<br>

## <a name="sim5polarization"></a>sim5polarization.c - Polarization properties of radiation



These routines help to calculate the change in polarization properties of radiation as it passes through a curved spacetime.

The orientation of the polarization vector (and thus the polarization angle) changes along the geodesic due to parallel transport. Knowing the inital direction of the vector, its orientation can be calculated at any place along the geodesic including the at-infinity limit. 


#### polarization_vector()

Photon polarization vector.


    DEVICEFUNC void polarization_vector(double k[4], sim5complex wp, sim5metric *metric, double f[4])

The returned polarization vector satisfies f.k=0 and f.f=1. Since f can be freely shifted by a multiple of k (f' = f + *k), it has a freedom in one compoment and it can always be set in such a way that its time-component is zero (f[0]=0). This routine returns f in such a form.

**Parameters**

* **k**: photon 4-momentum vector 
* **wp**: complex Walker-Penrose constant 
* **m**: metric 
* **f**: photon polarization vector (output)


**Return value**

Photon polarization vector f. 




<br>

---------

#### polarization_constant()

Walker-Penrose constant for a null geodesic.


    DEVICEFUNC sim5complex polarization_constant(double k[4], double f[4], sim5metric *metric)

Calculates components of Walker-Penrose (Walker&Penrose 1970) constant following Dexter(1916). returns Note: K1 and K2 relate to Connors,Piran,Stark(1980) kappa1=K2, kappa2=K1.

**Parameters**

* **k**: photon 4-momentum vector (satisfies k.k=0) 
* **f**: photon polarization vector (satisfies f.k=0) 
* **metric**: local metric


**Return value**

Complex Walker-Penrose constant $K_{wp} = K1 + I*K2$. 




<br>

---------

#### polarization_constant_infinity()

Walker-Penrose constant for a null geodesic.


    DEVICEFUNC sim5complex polarization_constant_infinity(double a, double alpha, double beta, double incl)

Calculates components of Walker-Penrose (Walker&Penrose 1970) constant following Dexter(1916). returns Note: K1 and K2 relate to Connors,Piran,Stark(1980) kappa1=K2, kappa2=K1.

**Parameters**

* **k**: photon 4-momentum vector (satisfies k.k=0) 
* **f**: photon polarization vector (satisfies f.k=0) 
* **metric**: local metric


**Return value**

Complex Walker-Penrose constant $K_{wp} = K1 + I*K2$. 




<br>

---------

#### polarization_angle_rotation()


    DEVICEFUNC double polarization_angle_rotation(double a, double inc, double alpha, double beta, sim5complex kappa)

<br>

---------

<br><br><br>

<br>

## <a name="sim5polyroots"></a>sim5polyroots.c - Polynomial roots



These routines give roots of quadratic, cubic and quartic equations calculated using analytical expressions. 


#### quadratic_eq()

Roots of qudratic equation.


    DEVICEFUNC int quadratic_eq(double pr, double pi, double qr, double qi, double zr[2], double zi[2])

Calculates the roots of $`z^2 + p z + q = 0`$ and the number of real roots, where p and q are complex numbers. The solution can have 2 real roots, 1 real and 1 complex roots, or 2 complex roots. The function uses the algorithm recommended in "Numerical Recipes".

**Parameters**

* **pr**: Re(p) 
* **pi**: Im(p) 
* **qr**: Re(q) 
* **qi**: Im(q) 
* **zr**: array to store real parts of quadratic roots 
* **zi**: array to store imaginary parts of quadratic roots


**Return value**

Number of real roots, plus the roots in zr and zi. 




<br>

---------

#### cubic_eq()

Roots of cubic equation.


    DEVICEFUNC int cubic_eq(double p, double q, double r, double zr[3], double zi[3])

Calculates the roots of $`z^3 + p z^2 + q z + r = 0`$ and the number of real roots, where p, q, r are real numbers.

**Parameters**

* **p**: coeficient p 
* **q**: coeficient q 
* **r**: coeficient r 
* **zr**: array to store real parts of cubic roots 
* **zi**: array to store imaginary parts of cubic roots


**Return value**

Number of real roots. The actual roots are returned in zr[] and zi[] arrays. 




<br>

---------

#### quartic_eq()

Roots of quartic equation.


    DEVICEFUNC int quartic_eq(double a3, double a2, double a1, double a0, double zr[4], double zi[4])

Calculates the roots of $` z^4+ a3 z^3 + a2 z^2 + a1 z + a0 = 0 `$ and the number of real roots, where a0, a1, a2, a3 are real numbers.

**Parameters**

* **a3**: coeficient a3 
* **a2**: coeficient a2 
* **a1**: coeficient a1 
* **a0**: coeficient a0 
* **zr**: array to store real parts of cubic roots 
* **zi**: array to store imaginary parts of cubic roots


**Return value**

Number of real roots, plus the roots in zr and zi. 




<br>

---------

#### quartic_eq_c()

Roots of quartic equation (variant).


    DEVICEFUNC void quartic_eq_c(double a3, double a2, double a1, double a0, int *nr, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4)

Calculates the roots of z^4+ a3 z^3 + a2 z^2 + a1 z + a0 = 0 and the number of real roots, where a0, a1, a2, a3 are real numbers. This variant returns all the root directly in 4 separate variables.

**Parameters**

* **a3**: coeficient a3 
* **a2**: coeficient a2 
* **a1**: coeficient a1 
* **a0**: coeficient a0 
* **nr**: number or real roots (output) 
* **z1**: first root (output) 
* **z2**: second root (output) 
* **z3**: third root (output) 
* **z4**: fourth root (output)


**Return value**

Number of real roots and the roots in nrr and z1...z4. 




<br>

---------

<br><br><br>

<br>

## <a name="sim5radiation"></a>sim5radiation.c - Radiative processes routines



Provides routines for radiative processes. 


#### blackbody_Iv()

Specific radiance of black-body radiation.


    DEVICEFUNC double blackbody_Iv(double T, double hardf, double cos_mu, double E)

Gives specific intensity {dE}{dt dA d dE}] at energy `E` of radiation of a black body of temperatute `T`. Assumes either isotropic or limb-darkened emission and a color correction of effective temperature by factor `hardf`.

**Parameters**

* **T**: temperature [K] 
* **hardf**: hardening factor (effective temperature correction) 
* **cos_mu**: cosine of emission direction with respect to the normal to the emission surface; set cos_mu>=0 for limb-darkened emission and cos_mu<0 for isotropic emission 
* **E**: energy [keV]


**Return value**

specific intensity in units [erg cm^-2 s^-1 keV^-1 srad^-1] 




<br>

---------

#### blackbody()

Specific radiance of black-body radiation.


    DEVICEFUNC void blackbody(double T, double hardf, double cos_mu, double E[], double Iv[], int N)

Gives specific intensity $`\frac{dE}{dt dA d\Omega dE}`$ at energy `E` of radiation of a black body of temperatute `T`. Assumes either isotropic or limb-darkened emission and a color correction of effective temperature by factor `hardf`.

**Parameters**

* **T**: temperature [K] 
* **hardf**: hardening factor (effective temperature correction) 
* **cos_mu**: cosine of emission direction with respect to the normal to the emission surface; set cos_mu>=0 for limb-darkened emission and cos_mu<0 for isotropic emission 
* **E**: array of energies (input) [keV] 
* **Iv**: array of specific intensities (output) [erg cm^-2 s^-1 keV^-1 srad^-1] 
* **N**: dimension of the E[] and Iv[] arrays


**Return value**

Iv[] array contains specific intensities for given energies 




<br>

---------

#### blackbody_photons()

Specific photon intensity of black-body radiation.


    DEVICEFUNC INLINE double blackbody_photons(double T, double hardf, double cos_mu, double E)

Same as `blackbody_Iv()` function, except it gives specific photon intensity $`\frac{dN}{dt dA d\Omega dE}`$.

**Return value**

specific intensity in units [photons cm^-2 s^-1 keV^-1 srad^-1] 




<br>

---------

#### blackbody_photons_total()

Number of photons that are emitted by a black-body surface.


    DEVICEFUNC double blackbody_photons_total(double T, double hardf)

Gives total number of photons $`\frac{dN}{dt dA d\Omega}`$ of all energies that are emitted form a unit surface of a black body of temperatute `T` into a unit solid angle. Assumes either isotropic or limb-darkened emission and a color correction of effective temperature by factor `hardf`.

**Parameters**

* **T**: temperature [K] 
* **hardf**: hardening factor (effective temperature correction) 
* **cos_mu**: cosine of emission direction with respect to the normal to the emission surface; set cos_mu>=0 for limb-darkened emission and cos_mu<0 for isotropic emission


**Return value**

number of photons [photons cm^-2 s^-1 srad^-1] 




<br>

---------

#### blackbody_photon_energy_random()

Draws a random photon energy that follows Planck distribution.


    DEVICEFUNC double blackbody_photon_energy_random(double T)

Picks a random photon energy of black-body radiation according to Planck energy distribution at given temperature. For derivation see http://arxiv.org/abs/1307.3635, part 3.3.1.

**Parameters**

* **T**: temperature of the distribution [K] (including hardening factor)


**Return value**

photon energy in [keV] 




<br>

---------

<br><br><br>

<br>

## <a name="sim5raytrace"></a>sim5raytrace.c - Raytracing



Provides routines for step-wise raytracing. 


| Name | Description | Value |
|------|-------------|-------|
| raytrace_max_error | default maximal relative error of raytracing | 1e-2 |
#### raytrace_prepare()

Raytracing with step-wise null-geodesic integration.


    DEVICEFUNC void raytrace_prepare(double bh_spin, double x[4], double k[4], double presision_factor, int options, raytrace_data *rtd)

Makes one step along the geodesic that is tangent to `k` and updates input vectors with new values. The integration method follows Dolence+2009 (http://adsabs.harvard.edu/abs/2009ApJS..184..387D).

The routine automatically controls size of the step based on curvature and required precission. On each call to raytrace() the routine takes a step size, which is smaller of the two: the internally chosen step size and step size that is passed on input in `step`.

Numerical precision of integration is driven by the precision_factor modifier and raytrace_max_error constant; roughly, final error (after whole geodesic is integrated) is: maximal relative error = (a factor of few) * raytrace_max_error * precision_factor.

**Parameters**

* **bh_spin**: black hole spin 
* **x**: initial position vector 
* **k**: initial direction vector (photon 4-momentum) 
* **f**: initial polarization vector (optional, can be NULL) 
* **precision_factor**: precision factor 
* **options**: additional options 
* **rtd**: raytracing data 




<br>

---------

#### raytrace_error()

Raytracing error.


    DEVICEFUNC double raytrace_error(double x[4], double k[4], raytrace_data *rtd)

Gives relative error in raytracing in terms of relative difference of Carter's constant. Useful for checking precission of integration.

**Parameters**

* **x**: position vector 
* **k**: direction vector 
* **f**: polarization vector (optional, can be NULL) 
* **rtd**: raytracing data


**Return value**

Relative error in raytracing. 




<br>

---------

<br><br><br>

<br>

## <a name="sim5roots"></a>sim5roots.c - Root finding



Routines for finding roots of functions numericaly. 


#### rtbis()

Root finding by bisection.


    long rtbis(double x1, double x2, double xacc, double(*fx)(double), double *result)

Finds root of a function on an interval. Using bisection method, it finds the root of a function `fx` that is known to lie between `x1` and `x2`. The root, returned as `result`, will be refined until its accuracy is +/- xacc.

**Parameters**

* **x1**: left boundary of the interval where the root is searched for 
* **x2**: right boundary of the interval where the root is searched for 
* **xacc**: accuracy 
* **fx**: function


**Return value**

Returns 1 if OK and the root position in result, 0 if error. 




<br>

---------

<br><br><br>

<br>

## <a name="sim5utils"></a>sim5utils.c - Utility routines




| Name | Description | Value |
|------|-------------|-------|
| MAXTOKENS |  | 1024 |
| MAXLINE |  | 8096 |
| MINLEN |  | 3 |
#### gprintf()


    void gprintf(FILE *file, const char *templatex,...)

<br>

---------

#### warning()


    void warning(const char *templatex,...)

<br>

---------

#### error()


    void error(const char *templatex,...)

<br>

---------

#### sort_array_d_compare_func()


    int sort_array_d_compare_func(const void *x, const void *y)

<br>

---------

#### sort_array()


    void sort_array(double *array, int N)

<br>

---------

#### sort_array_f_compare_func()


    int sort_array_f_compare_func(const void *x, const void *y)

<br>

---------

#### sort_array_f()


    void sort_array_f(float *array, int N)

<br>

---------

#### array_alloc()


    void * array_alloc(size_t capacity, size_t element_size)

<br>

---------

#### array_free()


    void array_free(void *ptr)

<br>

---------

#### array_realloc()


    void * array_realloc(void *array, size_t new_capacity)

<br>

---------

#### array_count()


    long array_count(void *array)

<br>

---------

#### array_capa()


    long array_capa(void *array)

<br>

---------

#### array_esize()


    size_t array_esize(void *array)

<br>

---------

#### array_push()


    void array_push(void **array_ptr, const void *data)

<br>

---------

#### array_push_int()


    void array_push_int(void **array_ptr, const int data)

<br>

---------

#### array_push_long()


    void array_push_long(void **array_ptr, const long data)

<br>

---------

#### array_push_float()


    void array_push_float(void **array_ptr, const float data)

<br>

---------

#### array_push_double()


    void array_push_double(void **array_ptr, const double data)

<br>

---------

#### array_exists()


    int array_exists(void *array, const void *data)

<br>

---------

#### array_push_if_not_exists()


    void array_push_if_not_exists(void **array_ptr, const void *data)

<br>

---------

#### array_reverse()


    void array_reverse(void *array)

<br>

---------

#### key_value_get()


    char * key_value_get(const char *string, const char *key)

<br>

---------

#### split()


    char ** split(char *string, char *delim)

<br>

---------

#### getlinecount()


    long getlinecount(FILE *f)

<br>

---------

#### backtrace()


    void backtrace()

<br>

---------

<br><br><br>

