# SIM5 Library Reference

<br/>

## Content
* [Constants and unit conversions (sim5const)](#sim5const)
* [Thin disk routines (sim5disk-nt)](#sim5disk-nt)
* [ (sim5distributions)](#sim5distributions)
* [ (sim5integration)](#sim5integration)
* [Interpolation routines (sim5interpolation)](#sim5interpolation)
* [ (sim5kerr-geod)](#sim5kerr-geod)
* [Basic spacetime routines (sim5kerr)](#sim5kerr)
* [ (sim5math)](#sim5math)
* [ (sim5polyroots)](#sim5polyroots)
* [Radiative processes routines (sim5radiation)](#sim5radiation)
* [Raytracing (sim5raytrace)](#sim5raytrace)
* [Root finding (sim5roots)](#sim5roots)
* [ (sim5utils)](#sim5utils)


<br/>

## <a name="sim5const"></a>sim5const.h - Constants and unit conversions



The library's default unit system is CGS. Sometimes constants in SI units (prefixed by `si_`) or in geometrical units (G=c=1; prefixed by `gu_`) are needed.

To convert between different units, some frequent conversion factors are defined. E.g. by multipying an energy value given in ergs by `erg2kev` factor makes a conversion to kiloelectronvolts. 


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
| atomig_mass_unit | atomic mass unit [g] | 1.660539e−24 |
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
<br/><br/><br/>

## <a name="sim5disk-nt"></a>sim5disk-nt.c - Thin disk routines



Provides routines for Novikov-Thorne thin disk model 


#### disk_nt_setup()

Set up a relativistic (Novikov-Thorne) model of thin disk.


    DEVICEFUNC int disk_nt_setup(double M, double a, double mdot_or_L, double alpha, int _options)

The disk can be set up using either mass accretion rate or luminosity. In case of Mdot, this is specified as a ratio of actual accretion rate in grams per second to the Eddington mass acretion rate corresponding to the given mass. The Eddington mass accretion rate constant is declared in sim5const.h. In case of luminosity, the model sets accretion rate such that the integrated disk luminosity matches the given value in ergs/sec relative to the Eddington luminosity again declared in sim5const.h.

**Parameters**

* **M**: mass of central BH [M_sun] 
* **a**: spin of central BH [0..1] 
* **mdot_or_L**: mass accretion rate (default) or luminosity (both in eddington units; see sim5const.h) 
* **alpha**: viscosity parameter 
* **options**: optional switches (default=0; switches can be combined with `+` operator)
  * DISK_NT_OPTION_LUMINOSITY: `mdot_or_L` parameter is interpreted as luminosity





**Return value**

A status code (currently zero) 




---------
<br/>

#### disk_nt_finish()

Finalize the disk model.


    DEVICEFUNC void disk_nt_finish()

Cleans up and frees memory. 

---------
<br/>

#### disk_nt_r_min()

Minimal radius of the disk (disk inner edge).


    DEVICEFUNC double disk_nt_r_min()

Provides minimal value for radius for which the functions provide valid results. For NT disk, this corresponds to the radius of marginally stable orbit (r_ms), where there is zero torque in the fluid.

**Return value**

Radius of disk inner edge [GM/c2] 




---------
<br/>

#### disk_nt_flux()

Local flux from one side of the disk.


    DEVICEFUNC double disk_nt_flux(double r)

Provides radial radiation flux dependence for Novikov-Thorne accretion disk. Formulae based on Page&Thorne(1974) http://adsabs.harvard.edu/abs/1974ApJ...191..499P

Note the retuned flux is local flux, i.e. flux measured by observer that is at rest with respect to the fluid.

**Parameters**

* **r**: radius of emission [GM/c2]


**Return value**

Total outgoing flux from unit area on one side of the disk [erg cm-2 s-1]. 




---------
<br/>

#### disk_nt_lum()

Gets total disk luminosity.


    DEVICEFUNC double disk_nt_lum()

Luminosity is obtained by integrating local flux over the surface area of the disk going into the whole sky (4pi solid angle). The integration makes a proper transformation of local flux to coordinate frame, but ignores other relativistic effects like light bending.



$$L = 2 * 2\pi \int F(r) r dr$$

**Return value**

Total disk luminosity of both surfaces [erg s-1] 




---------
<br/>

#### disk_nt_mdot()

Mass accretion rate.


    DEVICEFUNC double disk_nt_mdot()

Returns mass accretion rate in Eddington units. See `disk_nt_setup()` for details.

**Return value**

Mass accretion rate in Eddington units. 




---------
<br/>

#### disk_nt_temp()

Effective temperature.


    DEVICEFUNC double disk_nt_temp(double r)

Returns disk effective temperature at given radius. The temperature value corresponds to total outgoing flux at thar radius through the Steffan-Boltzmann law.

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Effective surface temperatute in [K]. 




---------
<br/>

#### disk_nt_sigma()

Column density.


    DEVICEFUNC double disk_nt_sigma(double r)

Returns midplane column density of the fluid at given radius for the first two zones according to formulae from Chandrasekhar book.



$$ \Sigma = \int_0^H \rho dz $$

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Midplane column density in [g/cm2]. 




---------
<br/>

#### disk_nt_ell()

Specific angular momentum.


    DEVICEFUNC double disk_nt_ell(double r)

Returns specific angular momentum of the fluid at given radius.

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Specific angular momentum in [g.u.]. 




---------
<br/>

#### disk_nt_vr()

Radial velocity.


    DEVICEFUNC double disk_nt_vr(double r)

Returns bulk radial velocity of the fluid at given radius. (For thin disk approximation, this is always zero.)

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Radial velocity in [speed_of_light]. 




---------
<br/>

#### disk_nt_h()

Surface height.


    DEVICEFUNC double disk_nt_h(double r)

Returns height of the surface of the disk (effective photosphere) above midplane at given radius. (For thin disk approximation, this is always zero.)

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Radial velocity in [speed_of_light]. 




---------
<br/>

#### disk_nt_dhdr()

Derivative of surface height.


    DEVICEFUNC double disk_nt_dhdr(double r)

Returns surface profile as derivative $dH/dR$ of its height above midplane at given radius. (For thin disk approximation, this is always zero.)

**Parameters**

* **r**: radius (measured in equatorial plane) [rg]


**Return value**

Derivative of surface height. 




---------
<br/>

#### disk_nt_dump()

Dump function printing disk provile.


    DEVICEFUNC void disk_nt_dump()

The function prints to stdout profile of all quantities as a function fo radius from r_ms to some outer radius (~2000 rg). 

---------
<br/>

#### disk_nt_find_mdot_for_luminosity()


    DEVICEFUNC double disk_nt_find_mdot_for_luminosity(double L0)

---------
<br/>

<br/><br/><br/>

## <a name="sim5distributions"></a>sim5distributions.c - 


#### distrib_init()

Creates a distribution based on given probability density function (PDF).


    DEVICEFUNC void distrib_init(sim5distrib *d, double(*pdf)(double), double x_min, double x_max, int N)

Uses given probability density function to initiate the internal data structure with calculated cummulative distribution its inverse.

**Parameters**

* **d**: pointer to structure that stores the disribution data 
* **pdf**: pointer to probability density function 
* **x_min**: left boundary for x 
* **x_min**: right boundary for x 
* **N**: number of samples to take over the interval [x_min:x_max] 




---------
<br/>

#### distrib_done()

Frees the internal data for the distribution.


    DEVICEFUNC void distrib_done(sim5distrib *d)

**Parameters**

* **d**: pointer to structure that stores the disribution data 




---------
<br/>

#### distrib_hit()

Generates a functional value on interval [x_min:x_max] according to the distribution.


    DEVICEFUNC INLINE double distrib_hit(sim5distrib *d)

With the help of precomputed cummulative distribution function it generates a value form the interval [x_min:x_max] and returns tnis value.

**Parameters**

* **d**: pointer to structure that stores the disribution data


**Return value**

value from the distribution 




---------
<br/>

<br/><br/><br/>

## <a name="sim5integration"></a>sim5integration.c - 


| Name | Description | Value |
|------|-------------|-------|
| NMAX |  | 23 |
| NMAX |  | 23 |
#### integrate_trapezoid_rule()

Integration core routine based on trapezoid rule.


    DEVICEFUNC void integrate_trapezoid_rule(double(*f)(double), double a, double b, int n, double *s)


* computes the `n`-th stage of refinement of an extended trapezoidal rule for function <f> and limits `a` and `b`
* when called with n=1, the routine returns `s` as the crudest estimate of $\int^a_b f(x) dx$
* subsequent calls with n=2,3,... (in that sequential order) will improve the accuracy of `s` by adding 2^(n-2) additional interior points.
* `s` must not be modified between sequential calls
* implementation is based on Numerical Recipes with the improvement by GSL, where a varible pointer is passed to the routine to store the value of the latest refinement instead of that being a global variable
* calling scheme: `for(j=1; j<=M+1; j++) trapezoid_rule(func, a, b, j, &answer);` 




---------
<br/>

#### integrate_trapezoid()

Integration of function using trapezoid rule.


    DEVICEFUNC double integrate_trapezoid(double(*f)(double), double a, double b, double acc)


* computes the integral ^a_b f(x) dx in a series of refinement steps until relative accuracy is better that <acc> or the maximum predefined number of steps is reached
* accuracy should not be increased beyond ~10^-6 as roundoff errors start to accumulate if too many steps are taken 




---------
<br/>

#### integrate_simpson()

Integration of function using simpson rule.


    DEVICEFUNC double integrate_simpson(double(*f)(double), double a, double b, double acc)


* computes the integral ^a_b f(x) dx in a series of refinement steps until relative accuracy is better that <acc> or the maximum predefined number of steps is reached
* simpson rule is generally more efficient than trapezoid rule when the function to be integrated has a finite 4th derivative (continuous 3rd derivative)
* accuracy should not be increased beyond ~10^-6 as roundoff errors start to accumulate if too many steps are taken 




---------
<br/>

#### gammln()


    double gammln(double xx)

---------
<br/>

#### gauleg()


    DEVICEFUNC void gauleg(double x1, double x2, double x[], double w[], int n)

---------
<br/>

#### qgaus()


    float qgaus(float(*func)(float), float a, float b)

---------
<br/>

#### qgaus_general()


    double qgaus_general(double w[], double y[], int N, double a, double b)

---------
<br/>

<br/><br/><br/>

## <a name="sim5interpolation"></a>sim5interpolation.c - Interpolation routines



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




---------
<br/>

#### sim5_interp_data_push()

Pushes data into interpolation object.


    DEVICEFUNC void sim5_interp_data_push(sim5interp *interp, double x, double y)

Adds a data point [x,y] into the interpolation object. This function is for filling interpolation object that has been created with option INTERP_DATA_BUILD with interolated data. The data must come in ordered sequence (x[i] < x[i+1])

**Parameters**

* **interp**: interpolation object 
* **x**: x-value of data point 
* **y**: y-value of data point 




---------
<br/>

#### sim5_interp_eval()

Interpolated data evaluation.


    DEVICEFUNC double sim5_interp_eval(sim5interp *interp, double x)

Makes the evalutaion on interpolated grid at given point.

**Parameters**

* **interp**: interpolation object 
* **x**: value for which to get interpolated value


**Return value**

Interpolated value. 




---------
<br/>

#### sim5_interp_done()

Interpolation finalization.


    DEVICEFUNC void sim5_interp_done(sim5interp *interp)

Frees the interpolation object interp (including a copied data, if necessary).

**Parameters**

* **interp**: interpolation object 




---------
<br/>

#### sim5_interp_alloc()

Alloc interpolation object memory.


    DEVICEFUNC sim5interp * sim5_interp_alloc()

Makes a memory allocation for interpolation object. This function should be used for heap-allocated variant of usage: sim5interp* interp;
interp = sim5_interp_alloc();
sim5_interp_init(interp, ...);
sim5_interp_free(interp);


**Return value**

Interpolation object. 




---------
<br/>

#### sim5_interp_free()

Free interpolation object memory.


    DEVICEFUNC void sim5_interp_free(sim5interp *interp)

Frees the interpolation object interp that had been previously alocated by `sim5_interp_alloc()`.

**Parameters**

* **interp**: interpolation object 




---------
<br/>

<br/><br/><br/>

## <a name="sim5kerr-geod"></a>sim5kerr-geod.c - 


| Name | Description | Value |
|------|-------------|-------|
| theta_int |  | (g->mK*jacobi_icn((x)/sqrt(g->m2p),g->m2)) |
| theta_inv |  | (sqrt(g->m2p)*jacobi_cn((x)/g->mK,g->m2)) |
#### geodesic_priv_R_roots()


    DEVICEFUNC int geodesic_priv_R_roots(geodesic *g, int *status)

---------
<br/>

#### geodesic_priv_T_roots()


    DEVICEFUNC int geodesic_priv_T_roots(geodesic *g, int *status)

---------
<br/>

#### geodesic_init_inf()

Makes setup for geodesic that is specified by impact parameters at infinity.


    DEVICEFUNC int geodesic_init_inf(double i, double a, double alpha, double beta, geodesic *g, int *status)

Parameters: i - inclination angle of observer (angle between BH rotation axis and direction to observer) [radians] a - BH spin [0..1] alpha - impact parameter in horizontal direction [GM/c^2] beta - impact parameter in vertical direction [GM/c^2] g - structure with information about geodesic status - status variable

Returns:
* status is 1 if setup is sucessfull, 0 in case of an error
* information about geodesic is stored in structure g 




---------
<br/>

#### geodesic_init_src()

Makes setup for geodesic that is specified by a point and direction (4-momentum vector).


    DEVICEFUNC int geodesic_init_src(double a, double r, double m, double k[4], int bpa, geodesic *g, int *status)

Parameters: a - BH spin [0..1] r - radial coordinate [GM/c^2] m - poloidal coordinate [cos(theta)] k - direction of the photon (4-momentum vector) bpa - position with respect to periastron (0=before periastron, 1=after periastron) g - structure with information about geodesic status - status variable

Returns:
* status is 1 if setup is sucessfull, 0 in case of an error
* information about geodesic is stored in structure g 




---------
<br/>

#### geodesic_P_int()

Returns the value of position integral along geodesics at point r.


    DEVICEFUNC double geodesic_P_int(geodesic *g, double r, int bpa)

The function gives the value of the integral P =  1/{R} dr =  1/{} d This integral is integrated from infinity to the point, where the geodesic reaches radius r either before or behind its periastron. The value of the integral increases monotonicly from infinity along the geodesic and thus we use it for parametrizing it. Note it is not affine parameter, which would be other choice for parametrization.

Parameters: g - geodesic r - radial coordinate [GM/c^2] bpa - position with respect to periastron (0=before periastron, 1=after periastron)

Returns: Value of the position integral between infinity and given point. 

---------
<br/>

#### geodesic_position()

Returns point on the geodesic, where the position integral gains value P.


    DEVICEFUNC void geodesic_position(geodesic *g, double P, double x[4])

The integral is evaluted in the form phi = a  (2r-al)/(Delta {R}) dr + l  sin^-2() / {Theta} d

Parameters: g - geodesic P - value of the position integral x[out] - coordinate

Returns: Fills x[] with position. 

---------
<br/>

#### geodesic_position_rad()

Gives radius at which the position integral gains value P.


    DEVICEFUNC double geodesic_position_rad(geodesic *g, double P)

Parameters: g - geodesic P - value of the position integral

Returns: Radius [GM/c^2] 

---------
<br/>

#### geodesic_position_pol()

Gives poloidal (theta) angle at which the position integral gains value P.


    DEVICEFUNC double geodesic_position_pol(geodesic *g, double P)

Parameters: g - geodesic P - value of the position integral

Returns: Cosine of theta angle. 

---------
<br/>

#### geodesic_position_azm()

Gives azimuthal (phi) angle at which the position integral gains value P.


    DEVICEFUNC double geodesic_position_azm(geodesic *g, double r, double m, double P)

The value of azimuthal angle is assumed to be zero at infinity and the function gives the change of the angle between the point [r,m] and infinity.

Parameters: g - geodesic P - value of the position integral

Returns: Phi angle in radians (can be more than 2pi) 

---------
<br/>

#### geodesic_timedelay()

Gives travel-time (timedelay) between positions P1 and P2.


    DEVICEFUNC double geodesic_timedelay(geodesic *g, double P1, double P2)

Parameters: g - geodesic P1 - value of the position integral at point A P2 - value of the position integral at point B

Returns: timedelay between position P1 and P1 (point A and B) 

---------
<br/>

#### geodesic_find_midplane_crossing()


    DEVICEFUNC double geodesic_find_midplane_crossing(geodesic *g, int order)

---------
<br/>

#### geodesic_follow()


    DEVICEFUNC void geodesic_follow(geodesic *g, double step, double *P, double *r, double *m, int *status)

---------
<br/>

#### geodesic_s2i_int_R()

Returns the value of position integral along geodesics at point r.


    DEVICEFUNC double geodesic_s2i_int_R(geodesic *g, double r, int bpa)

The function gives the value of the integral P =  1/{R} dr =  1/{} d This integral is integrated from infinity to the point, where the geodesic reaches radius r either before or behind its periastron. The value of the integral increases monotonicly from infinity along the geodesic and thus we use it for parametrizing it. Note it is not affine parameter, which would be other choice for parametrization.

Parameters: g - geodesic r - radial coordinate [GM/c^2] bpa - position with respect to periastron (0=before periastron, 1=after periastron)

Returns: Value of the position integral between infinity and given point. 

---------
<br/>

#### geodesic_s2i_inv_R()

Returns the value of position integral along geodesics at point r.


    DEVICEFUNC double geodesic_s2i_inv_R(geodesic *g, double r, int *bpa)

The function gives the value of the integral P =  1/{R} dr =  1/{} d This integral is integrated from infinity to the point, where the geodesic reaches radius r either before or behind its periastron. The value of the integral increases monotonicly from infinity along the geodesic and thus we use it for parametrizing it. Note it is not affine parameter, which would be other choice for parametrization.

Parameters: g - geodesic r - radial coordinate [GM/c^2] bpa - position with respect to periastron (0=before periastron, 1=after periastron)

Returns: Value of the position integral between infinity and given point. 

---------
<br/>

#### geodesic_s2i_int_T_eqplane()


    DEVICEFUNC double geodesic_s2i_int_T_eqplane(geodesic *g, int order, double *dk2dx)

---------
<br/>

#### geodesic_s2i_inv_T()


    DEVICEFUNC double geodesic_s2i_inv_T(geodesic *g, double T, double *dk2dx)

---------
<br/>

#### geodesic_s2i_timedelay()


    DEVICEFUNC double geodesic_s2i_timedelay(geodesic *g, double x, double *opt_r, double *opt_m)

---------
<br/>

#### geodesic_s2i_phi()


    DEVICEFUNC double geodesic_s2i_phi(geodesic *g, double x, double *opt_r, double *opt_m)

---------
<br/>

#### geodesic_s2i_afp()


    DEVICEFUNC double geodesic_s2i_afp(geodesic *g, double x, double *opt_r, double *opt_m)

---------
<br/>

#### geodesic_s2i_solution_eqplane()


    DEVICEFUNC void geodesic_s2i_solution_eqplane(geodesic *g, int order, double *r, int *beyond_pa, double *phi, double *time_delay, int *status)

---------
<br/>

<br/><br/><br/>

## <a name="sim5kerr"></a>sim5kerr.c - Basic spacetime routines



Provides basic routines for doing physics in Kerr and Minkowski spacetimes (metric, connection, tetrads, vector algebra, orbital motion, photon propagation, etc). 


#### flat_metric()

Flat (Minkowski) spacetime metric.


    DEVICEFUNC void flat_metric(double r, double m, sim5metric *metric)

Returns covariant Minkowski metric $\eta_\mu\nu$ in spherical coordinates (t,r,theta,phi).

**Parameters**

* **r**: radial coordinate 
* **m**: poloidal coordinate $m=\cos\theta$


**Return value**

Metric components are returned in metric parameter. 




---------
<br/>

#### flat_metric_contravariant()

Flat (Minkowski) spacetime metric (contravariant).


    DEVICEFUNC void flat_metric_contravariant(double r, double m, sim5metric *metric)

Returns contravariant Minkowski metric $\eta^\mu\nu$ in spherical coordinates (t,r,theta,phi).

**Parameters**

* **r**: radial coordinate 
* **m**: poloidal coordinate $m=\cos\theta$


**Return value**

Metric components are returned in metric parameter. 




---------
<br/>

#### kerr_metric()

Kerr spacetime metric.


    DEVICEFUNC void kerr_metric(double a, double r, double m, sim5metric *metric)

Returns Kerr metric $g_\mu\nu$ in spherical coordinates (t,r,theta,phi).

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate 
* **m**: poloidal coordinate $m=\cos\theta$


**Return value**

Metric components are returned in metric parameter. 




---------
<br/>

#### kerr_metric_contravariant()

Kerr spacetime metric (contravariant).


    DEVICEFUNC void kerr_metric_contravariant(double a, double r, double m, sim5metric *metric)

Returns contravariant Kerr metric $g^\mu\nu$ in spherical coordinates (t,r,theta,phi).

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate 
* **m**: poloidal coordinate $m=\cos\theta$


**Return value**

Metric components are returned in metric parameter. 




---------
<br/>

#### kerr_connection()

Christoffel symbol components for Kerr metric ( $Gamma^\mu_\alpha\beta$).


    DEVICEFUNC void kerr_connection(double a, double r, double m, double G[4][4][4])

Returns a matrix of connection coefficients for Kerr metric.

(!) NOTE: Christoffel tensor Gamma^i_(jk) is symmetric in lower two indices (j,k). This function only evaluates components of the tensor, where j<k and multiplies the value of these coeficients by a factor of two. In this way, both summation over j=1..4,k=1..4 and over j=1..4,k=j..4 will give the same result; however, care must be taken when summing Gamma^i_(jk) U^j V^k where factor 0.5 has to be used: `0.5*G[i][j][k]*(u[j]*v[k] + u[k]*v[j])`

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate 
* **m**: poloidal coordinate $m=\cos\theta$


**Return value**

Connection coeficients are returned in G parameter. 




---------
<br/>

#### flat_connection()

Christoffel symbol components for Minkowski metric ( $Gamma^\mu_\alpha\beta$).


    DEVICEFUNC void flat_connection(double r, double m, double G[4][4][4])

Returns a matrix of connection coefficients for Minkowski metric, i.e. Kerr metric in the limit of M=0 and a=0.

(!) NOTE: Christoffel tensor Gamma^i_(jk) is symmetric in lower two indices (j,k). This function only evaluates components of the tensor, where j<k and multiplies the value of these coeficients by a factor of two. In this way, both summation over j=1..4,k=1..4 and over j=1..4,k=j..4 will give the same result; however, care must be taken when summing Gamma^i_(jk) U^j V^k where factor 0.5 has to be used: `0.5*G[i][j][k]*(u[j]*v[k] + u[k]*v[j])`

**Parameters**

* **r**: radial coordinate 
* **m**: poloidal coordinate $m=\cos\theta$


**Return value**

Connection coeficients are returned in G parameter. 




---------
<br/>

#### Gamma()

Change of vector along a trajectory.


    DEVICEFUNC INLINE void Gamma(double G[4][4][4], double U[4], double V[4], double result[4])

Returns product of summation `-G^i_(jk) U^j V^k`. This is useful for calculating parallel transport of vector `V` along a trajectory specified by tangent vector `U`, as it gives the derivative $dV/d\lambda$. (Note the definition of G tensor components in symmetric indices.)

**Parameters**

* **G**: connection coefficients 
* **U**: tangent vector 
* **V**: transported vector


**Return value**

Change of transported vector $dV/d\lambda$. 




---------
<br/>

#### vector()

Make a 4-vector.


    DEVICEFUNC INLINE void vector(double x[4], double x0, double x1, double x2, double x3)

Returns a 4-vector that has given components (contravarient form $X^\mu$).

**Parameters**

* **x0-x3**: components of the vector


**Return value**

Vector is returned in x parameter. 




---------
<br/>

#### vector_covariant()

Covariant version of a vector.


    DEVICEFUNC INLINE void vector_covariant(double V1[4], double V2[4], sim5metric *m)

Converts a standard (contravariant) form of a vector to covariant form ( $X^\mu -> X_\mu$). Metric `m` can be NULL in which case flat metric is assumed.

**Parameters**

* **V1**: contravariant vector 
* **V2**: covariant vector (output) 
* **m**: metric


**Return value**

Transformed vector is returned in V2 parameter. 




---------
<br/>

#### vector_norm()

Norm of a vector.


    DEVICEFUNC INLINE double vector_norm(double V[4], sim5metric *m)

Returns norm of a vector V, i.e. sqrt(V*V). Metric `m` can be NULL in which case flat metric is assumed. NOTE: only works for space-like vectors, where V*V>0.

**Parameters**

* **V**: vector 
* **m**: metric


**Return value**

Norm of a vector. 




---------
<br/>

#### vector_3norm()

3-norm of a vector (in flat spacetime).


    DEVICEFUNC INLINE double vector_3norm(double V[4])

Returns a 3-norm of a vector V (includes only 'spatial' components) in flat (Minkowski) spacetime, i.e. sqrt( V[i]*V[i]), where i=1..3.

**Parameters**

* **V**: vector 
* **m**: metric


**Return value**

3-norm of a vector. 




---------
<br/>

#### vector_multiply()

Multiply a vector.


    DEVICEFUNC INLINE void vector_multiply(double V[4], double factor)

Multiplies the components of a vector by given factor.

**Parameters**

* **V**: vector 
* **factor**: multiplication factor


**Return value**

Modified vector V. 




---------
<br/>

#### vector_norm_to()

Normalizes a vector to given size.


    DEVICEFUNC INLINE void vector_norm_to(double V[4], double norm, sim5metric *m)

Changes the vector components such that the norm of the vector becomes `|norm|`, i.e. $V.V = \rm{norm}$. The vector `V` must be space-like if norm>0, and it must be time-like if norm<0 (no checks are performed to enfoce that). Metric `m` can be NULL in which case flat metric is assumed. NOTE: vector cannot be null vector (V*V=0), nor norm can be zero; use vector_norm_to_null() in that case.

**Parameters**

* **V**: vector to be normalized 
* **norm**: normalization the vector should have 
* **m**: metric


**Return value**

Changes the components of vector V. 




---------
<br/>

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




---------
<br/>

#### dotprod()

Scalar product.


    DEVICEFUNC INLINE double dotprod(double V1[4], double V2[4], sim5metric *m)

Calculates scalar product of two 4-vectors ( $U^\mu V^\nu g_{\mu\nu}$). NOTE: if metric is NULL, flat space metric is used

**Parameters**

* **V1**: first vector 
* **V1**: second vector 
* **m**: metric


**Return value**

The value of scalar product 




---------
<br/>

#### tetrad_zamo()

Tetrad of a zero angular momentum observer (ZAMO).


    DEVICEFUNC void tetrad_zamo(sim5metric *m, sim5tetrad *t)

Returns basis vectors for ZAMO (locally non-rotating) observer tetrad e_{(a)}^{}. Note the orientation of theta-vector that points in opposite direction than d/d(theta), i.e. at the equatorial plane the theta-vector points in upward direction.

**Parameters**

* **m**: metric 
* **t**: tetrad


**Return value**

Returns tetrad vectors in t. 




---------
<br/>

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




---------
<br/>

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




---------
<br/>

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




---------
<br/>

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




---------
<br/>

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




---------
<br/>

#### r_bh()

Black hole event horizon radius.


    DEVICEFUNC INLINE double r_bh(double a)

**Parameters**

* **a**: black hole spin


**Return value**

Event horizon radius in [rg]. 




---------
<br/>

#### r_ms()

Radius of marginally stable orbit.


    DEVICEFUNC INLINE double r_ms(double a)

**Parameters**

* **a**: black hole spin


**Return value**

Marginally stable orbit radius in [rg]. 




---------
<br/>

#### r_mb()

Radius of marginally bound orbit.


    DEVICEFUNC INLINE double r_mb(double a)

Marginally bound orbit is minimal radius of bound (E<0) circular and parabolic orbits. http://adsabs.harvard.edu/abs/1972ApJ...178..347B, eq. 2.19.

**Parameters**

* **a**: black hole spin


**Return value**

Marginally bound orbit radius in [rg]. 




---------
<br/>

#### r_ph()

Radius of photon orbit.


    DEVICEFUNC INLINE double r_ph(double a)

Photon orbit is radius of unstable circular photon orbit. http://adsabs.harvard.edu/abs/1972ApJ...178..347B, eq. 2.18.

**Parameters**

* **a**: black hole spin


**Return value**

Marginally bound orbit radius in [rg]. 




---------
<br/>

#### OmegaK()

Angular frequency of Keplerian orbital motion.


    DEVICEFUNC INLINE double OmegaK(double r, double a)

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin


**Return value**

Angular frequency [Hz]. 




---------
<br/>

#### ellK()

Specific angular momentum of Keplerian orbital motion.


    DEVICEFUNC INLINE double ellK(double r, double a)

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin


**Return value**

Keplerian specific angular momentum [g.u.]. 




---------
<br/>

#### omega_r()

Angular frequency of radial epyciclic motion.


    DEVICEFUNC INLINE double omega_r(double r, double a)

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin


**Return value**

Angular velocity [Hz]. 




---------
<br/>

#### omega_z()

Angular frequency of vertical epyciclic motion.


    DEVICEFUNC INLINE double omega_z(double r, double a)

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin


**Return value**

Angular velocity [Hz]. 




---------
<br/>

#### Omega_from_ell()

Angular frequency corresponding to given angular momentum.


    DEVICEFUNC INLINE double Omega_from_ell(double ell, sim5metric *m)

**Parameters**

* **ell**: specific angular momentum [g.u.] 
* **m**: metric


**Return value**

Angular frequency [Hz]. 




---------
<br/>

#### ell_from_Omega()

Specific angular momentum corresponding to given angular frequency.


    DEVICEFUNC INLINE double ell_from_Omega(double Omega, sim5metric *m)

**Parameters**

* **Omega**: angular frequency [Hz] 
* **m**: metric


**Return value**

Specific angular momentum [g.u.]. 




---------
<br/>

#### gfactorK()

Redshift factor.


    double gfactorK(double r, double a, double l)

Relativistic correction to energy of photon that is emitted by fluid at Keplerian rotation in equatorial plane; includes Doppler effect and gravitational redshift.

**Parameters**

* **r**: radius [rg] 
* **a**: black hole spin 
* **l**: photon motion constant lambda


**Return value**

Redshift factor. 




---------
<br/>

#### photon_momentum()

Photon 4-momentum vector .


    DEVICEFUNC void photon_momentum(double a, double r, double m, double l, double q2, double r_sign, double m_sign, double k[4])

Returns photon 4-momentum vector k^ such that k*k=0 (null vector). The orientation of the resulting vector `k` is given by the signs of `r_sign` and `m_sign` parameters. See MTW; Rauch & Blandford (1994), Appendix A; Li+05

**Parameters**

* **a**: black hole spin 
* **r**: radial coordinate [rg] 
* **m**: poloidal coordinate [cos(theta)] 
* **l**: photon motion constant lambda 
* **q2**: photon motion constant Q^2 (Carter's constant) 
* **r_sign**: sign of k[1] component of resulting momentum vector 
* **m_sign**: sign of k[2] component of resulting momentum vector 
* **k**: resulting momentum vector (output)


**Return value**

Photon momentum vector k. 




---------
<br/>

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




---------
<br/>

#### photon_carter_const()

Carter's constant Q for null geodesic.


    DEVICEFUNC double photon_carter_const(double k[4], sim5metric *metric)

**Parameters**

* **k**: photon 4-momentum vector 
* **m**: metric


**Return value**

Carter's constants Q. 




---------
<br/>

#### photon_wp_const()

Walker-Penrose constant of null geodesic.


    DEVICEFUNC sim5complex photon_wp_const(double k[4], double f[4], sim5metric *metric)

Calculates components of Walker-Penrose (Walker&Penrose 1970) constant following Connors, Piran, Stark (1980): kappa_wp = kappa1 + I*kappa2 = (A1 - I*A2)*(r - I*a*cos(theta)) => kappa1 = +r*A1 - a*cos(theta)*A2; kappa2 = -r*A2 - a*cos(theta)*A1 returns kappa_wp = kappa1 + I*kappa2 Note the definition of kappa1 & kappa2, which is opposite to CPS(1980).

**Parameters**

* **k**: photon 4-momentum vector (satisfies k.k=0) 
* **f**: photon polarization vector (satisfies f.k=0) 
* **m**: metric


**Return value**

Complex Walker-Penrose constant. 




---------
<br/>

#### photon_polarization_vector()

Photon polarization vector.


    DEVICEFUNC void photon_polarization_vector(double k[4], sim5complex wp, sim5metric *metric, double f[4])

The returned polarization vector satisfies f.k=0 and f.f=0. Since f can be freely shifted by a multiple of k (f' = f + *k), it has a freedom in one compoment and it can always be set in such a way that its time-component is zero (f[0]=0). This routine returns f in such a form.

**Parameters**

* **k**: photon 4-momentum vector 
* **wp**: complex Walker-Penrose constant 
* **m**: metric 
* **f**: photon polarization vector (output)


**Return value**

Photon polarization vector f. 




---------
<br/>

#### fourvelocity_zamo()

4-velocity of ZAMO (locally non-rotating) observer.


    DEVICEFUNC INLINE void fourvelocity_zamo(sim5metric *m, double U[4])

**Parameters**

* **m**: metric 
* **U**: 4-velocity (output)


**Return value**

4-velocity U. 




---------
<br/>

#### fourvelocity_azimuthal()

4-velocity of azimuthally rotating observer.


    DEVICEFUNC INLINE void fourvelocity_azimuthal(double Omega, sim5metric *m, double U[4])

**Parameters**

* **Omega**: angular frequency of circular motion 
* **m**: metric 
* **U**: 4-velocity (output)


**Return value**

4-velocity U. 




---------
<br/>

#### fourvelocity_radial()

4-velocity of radially moving observer.


    DEVICEFUNC INLINE void fourvelocity_radial(double vr, sim5metric *m, double U[4])

**Parameters**

* **vr**: radial component of 4-velocity (dr/dtau) 
* **m**: metric 
* **U**: 4-velocity (output)


**Return value**

4-velocity U. 




---------
<br/>

#### ortho_tetrad_U()


    DEVICEFUNC void ortho_tetrad_U(double U[4], double g00, double g11, double g22, double g33, double g03, double e0[4], double e1[4], double e2[4], double e3[4])

---------
<br/>

#### ortho_tetrad_U_phi_r_motion()


    DEVICEFUNC void ortho_tetrad_U_phi_r_motion(double U[4], double g00, double g11, double g22, double g33, double g03, double e0[4], double e1[4], double e2[4], double e3[4])

---------
<br/>

#### fourvelocity_norm()


    DEVICEFUNC INLINE double fourvelocity_norm(double U1, double U2, double U3, double g00, double g11, double g22, double g33, double g03)

---------
<br/>

#### kappa_pw()


    DEVICEFUNC void kappa_pw(double a, double r, double m, double k[4], double f[4], double *kappa1, double *kappa2)

---------
<br/>

#### stokes_infty()


    DEVICEFUNC void stokes_infty(double a, double inc, double alpha, double beta, double kappa1, double kappa2, double *pol_angle)

---------
<br/>

<br/><br/><br/>

## <a name="sim5math"></a>sim5math.c - 


#### sim5round()


    DEVICEFUNC INLINE long sim5round(double num)

---------
<br/>

#### factorial()


    DEVICEFUNC INLINE long int factorial(long int n)

---------
<br/>

#### reduce_angle_pi()


    DEVICEFUNC INLINE double reduce_angle_pi(double phi)

---------
<br/>

#### reduce_angle_2pi()


    DEVICEFUNC INLINE double reduce_angle_2pi(double phi)

---------
<br/>

#### ensure_range()


    DEVICEFUNC int ensure_range(double *val, double min, double max, double acc)

---------
<br/>

#### sim5seed()


    DEVICEFUNC INLINE void sim5seed()

---------
<br/>

#### sim5rand()


    DEVICEFUNC INLINE unsigned long long sim5rand()

---------
<br/>

#### sim5urand()


    DEVICEFUNC INLINE double sim5urand()

---------
<br/>

#### cartesian2spherical1()


    DEVICEFUNC void cartesian2spherical1(double x, double y, double z, double Vx, double Vy, double Vz, double *Vr, double *Vh, double *Vf)

---------
<br/>

#### cartesian2spherical2()


    DEVICEFUNC void cartesian2spherical2(double cos_h, double sin_f, double cos_f, double Vx, double Vy, double Vz, double *Vr, double *Vh, double *Vf)

---------
<br/>

#### makeComplex()


    DEVICEFUNC INLINE sim5complex makeComplex(double r, double i)

---------
<br/>

#### nullComplex()


    DEVICEFUNC INLINE sim5complex nullComplex()

---------
<br/>

<br/><br/><br/>

## <a name="sim5polyroots"></a>sim5polyroots.c - 


#### quadratic_eq()

Roots of qudratic equation.


    DEVICEFUNC int quadratic_eq(double pr, double pi, double qr, double qi, double zr[2], double zi[2])

Calculates the roots of z^2 + p z + q = 0 and the number of real roots, where p and q are complex numbers. The solution can have 2 real roots, 1 real and 1 complex roots, or 2 complex roots. The function uses the algorithm recommended in "Numerical Recipes".

**Parameters**

* **pr**: Re(p) 
* **pi**: Im(p) 
* **qr**: Re(q) 
* **qi**: Im(q) 
* **zr**: array to store real parts of quadratic roots 
* **zr**: array to store imaginary parts of quadratic roots


**Return value**

Number of real roots, plus the roots in zr and zi. 




---------
<br/>

#### cubic_eq()

Roots of cubic equation.


    DEVICEFUNC int cubic_eq(double p, double q, double r, double zr[3], double zi[3])

Calculates the roots of z^3 + p z^2 + q z + r = 0 and the number of real roots, where p, q, r are real numbers.

**Parameters**

* **p**: coeficient p 
* **q**: coeficient q 
* **r**: coeficient r 
* **zr**: array to store real parts of cubic roots 
* **zr**: array to store imaginary parts of cubic roots


**Return value**

Number of real roots, plus the roots in zr and zi. 




---------
<br/>

#### sort_roots_re()


    DEVICEFUNC void sort_roots_re(double *r1, double *r2, double *r3, double *r4)

---------
<br/>

#### sort_mix()


    DEVICEFUNC void sort_mix(double *r1, double *r2, int *s)

---------
<br/>

#### sort_mix2()


    DEVICEFUNC void sort_mix2(double *r1, double *r2, int *s)

---------
<br/>

#### sort_roots()


    DEVICEFUNC void sort_roots(int *s, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4)

---------
<br/>

#### quartic_eq()

Roots of quartic equation.


    DEVICEFUNC int quartic_eq(double a3, double a2, double a1, double a0, double zr[4], double zi[4])

Calculates the roots of z^4+ a3 z^3 + a2 z^2 + a1 z + a0 = 0 and the number of real roots, where a0, a1, a2, a3 are real numbers.

**Parameters**

* **a3**: coeficient a3 
* **a2**: coeficient a2 
* **a1**: coeficient a1 
* **a0**: coeficient a0 
* **zr**: array to store real parts of cubic roots 
* **zr**: array to store imaginary parts of cubic roots


**Return value**

Number of real roots, plus the roots in zr and zi. 




---------
<br/>

#### quartic_eq_c()

Roots of quartic equation (variant).


    DEVICEFUNC void quartic_eq_c(double a3, double a2, double a1, double a0, int *nr, sim5complex *z1, sim5complex *z2, sim5complex *z3, sim5complex *z4)

Calculates the roots of z^4+ a3 z^3 + a2 z^2 + a1 z + a0 = 0 and the number of real roots, where a0, a1, a2, a3 are real numbers. This variant returns all the root directly in 4 separate variables.

**Parameters**

* **a3**: coeficient a3 
* **a2**: coeficient a2 
* **a1**: coeficient a1 
* **a0**: coeficient a0 
* **nrr**: number or real roots (output) 
* **z1**: first root (output) 
* **z2**: second root (output) 
* **z3**: third root (output) 
* **z4**: fourth root (output)


**Return value**

Number of real roots and the roots in nrr and z1...z4. 




---------
<br/>

<br/><br/><br/>

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




---------
<br/>

#### blackbody()

Specific radiance of black-body radiation.


    void blackbody(double T, double hardf, double cos_mu, double E[], double Iv[], int N)

Gives specific intensity $\frac{dE}{dt dA d\Omega dE}$ at energy `E` of radiation of a black body of temperatute `T`. Assumes either isotropic or limb-darkened emission and a color correction of effective temperature by factor `hardf`.

**Parameters**

* **T**: temperature [K] 
* **hardf**: hardening factor (effective temperature correction) 
* **cos_mu**: cosine of emission direction with respect to the normal to the emission surface; set cos_mu>=0 for limb-darkened emission and cos_mu<0 for isotropic emission 
* **E**: array of energies (input) [keV] 
* **Iv**: array of specific intensities (output) [erg cm^-2 s^-1 keV^-1 srad^-1] 
* **N**: dimension of the E[] and Iv[] arrays


**Return value**

Iv[] array contains specific intensities for given energies 




---------
<br/>

#### blackbody_photons()

Specific photon intensity of black-body radiation.


    DEVICEFUNC INLINE double blackbody_photons(double T, double hardf, double cos_mu, double E)

Same as `blackbody_Iv()` function, except it gives specific photon intensity $\frac{dN}{dt dA d\Omega dE}$.

**Return value**

specific intensity in units [photons cm^-2 s^-1 keV^-1 srad^-1] 




---------
<br/>

#### blackbody_photons_total()

Number of photons that are emitted by a black-body surface.


    DEVICEFUNC double blackbody_photons_total(double T, double hardf, double cos_mu)

Gives total number of photons $\frac{dN}{dt dA d\Omega}$ of all energies that are emitted form a unit surface of a black body of temperatute `T` into a unit solid angle. Assumes either isotropic or limb-darkened emission and a color correction of effective temperature by factor `hardf`.

**Parameters**

* **T**: temperature [K] 
* **hardf**: hardening factor (effective temperature correction) 
* **cos_mu**: cosine of emission direction with respect to the normal to the emission surface; set cos_mu>=0 for limb-darkened emission and cos_mu<0 for isotropic emission


**Return value**

number of photons [photons cm^-2 s^-1 srad^-1] 




---------
<br/>

#### blackbody_photon_energy_random()

Draws a random photon energy that follows Planck distribution.


    DEVICEFUNC double blackbody_photon_energy_random(double T)

Picks a random photon energy of black-body radiation according to Planck energy distribution at given temperature. For derivation see http://arxiv.org/abs/1307.3635, part 3.3.1.

**Parameters**

* **T**: temperature of the distribution [K] (including hardening factor)


**Return value**

photon energy in [keV] 




---------
<br/>

<br/><br/><br/>

## <a name="sim5raytrace"></a>sim5raytrace.c - Raytracing



Provides routines for step-wise raytracing. 


| Name | Description | Value |
|------|-------------|-------|
| raytrace_max_error | default maximal relative error of raytracing | 1e-2 |
#### raytrace_prepare()

Raytracing with step-wise null-geodesic integration.


    DEVICEFUNC void raytrace_prepare(double bh_spin, double x[4], double k[4], double f[4], double presision_factor, int options, raytrace_data *rtd)

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




---------
<br/>

#### raytrace()

Raytracing with step-wise null-geodesic integration.


    DEVICEFUNC void raytrace(double x[4], double k[4], double f[4], double *step, raytrace_data *rtd)

Makes one step along the geodesic that is tangent to `k` and updates input vectors with new values. The integration method follows Dolence+2009 (http://adsabs.harvard.edu/abs/2009ApJS..184..387D).

The routine automatically controls size of the step based on curvature and required precission. On each call to raytrace() the routine takes a step size, which is smaller of the two: the internally chosen step size and step size that is passed on input in `step`.

Numerical precision is driven by the precision_factor modifier [see raytrace_prepare()]; rtd->error gives the error in the current step and it should be bellow raytrace_max_error*1e-2; so it is adviceable to check rtd->error continuously after each step and stop integration when the error goes above ~1e-3; at the end of integration one should then check relative difference in Carter's constant with `raytrace_error()`.

**Parameters**

* **x**: position vector 
* **k**: direction vector (photon 4-momentum) 
* **f**: polarization vector (optional, can be NULL) 
* **step**: on input gives maximal step that the routine can take [GM/c2]; on output gives size of step that has actually been taken 
* **rtd**: raytracing data


**Return value**

Position, direction and polarization (if not null) vectors and step size are updated, rtd has the relative error. 




---------
<br/>

#### raytrace_error()

Raytracing error.


    DEVICEFUNC double raytrace_error(double x[4], double k[4], double f[4], raytrace_data *rtd)

Gives relative error in raytracing in terms of relative difference of Carter's constant. Useful for checking precission of integration.

**Parameters**

* **x**: position vector 
* **k**: direction vector 
* **f**: polarization vector (optional, can be NULL) 
* **rtd**: raytracing data


**Return value**

Relative error in raytracing. 




---------
<br/>

<br/><br/><br/>

## <a name="sim5roots"></a>sim5roots.c - Root finding



Provides routines for ffinding roots of functions. 


#### rtbis()

Root finding.


    long rtbis(double x1, double x2, double xacc, double(*fx)(double), double *result)

Finds root of a function on an interval. Using bisection method, finds the root of a function `fx` that is known to lie between `x1` and `x2`. The root, returned as `result`, will be refined until its accuracy is +/- xacc.

**Parameters**

* **x1**: left boundary of the interval where the root is searched for 
* **x2**: right boundary of the interval where the root is searched for 
* **xacc**: accuracy 
* **fx**: function


**Return value**

Returns 1 if OK and the root position in result, 0 if error. 




---------
<br/>

<br/><br/><br/>

## <a name="sim5utils"></a>sim5utils.c - 


| Name | Description | Value |
|------|-------------|-------|
| MAXTOKENS |  | 1024 |
| MAXLINE |  | 8096 |
| MINLEN |  | 3 |
#### gprintf()


    void gprintf(FILE *file, const char *templatex,...)

---------
<br/>

#### warning()


    void warning(const char *templatex,...)

---------
<br/>

#### error()


    void error(const char *templatex,...)

---------
<br/>

#### sort_array_d_compare_func()


    int sort_array_d_compare_func(const void *x, const void *y)

---------
<br/>

#### sort_array()


    void sort_array(double *array, int N)

---------
<br/>

#### sort_array_f_compare_func()


    int sort_array_f_compare_func(const void *x, const void *y)

---------
<br/>

#### sort_array_f()


    void sort_array_f(float *array, int N)

---------
<br/>

#### array_alloc()


    void * array_alloc(size_t capacity, size_t element_size)

---------
<br/>

#### array_free()


    void array_free(void *ptr)

---------
<br/>

#### array_realloc()


    void * array_realloc(void *array, size_t new_capacity)

---------
<br/>

#### array_count()


    long array_count(void *array)

---------
<br/>

#### array_capa()


    long array_capa(void *array)

---------
<br/>

#### array_esize()


    size_t array_esize(void *array)

---------
<br/>

#### array_push()


    void array_push(void **array_ptr, const void *data)

---------
<br/>

#### array_push_int()


    void array_push_int(void **array_ptr, const int data)

---------
<br/>

#### array_push_long()


    void array_push_long(void **array_ptr, const long data)

---------
<br/>

#### array_push_double()


    void array_push_double(void **array_ptr, const double data)

---------
<br/>

#### array_exists()


    int array_exists(void *array, const void *data)

---------
<br/>

#### array_push_if_not_exists()


    void array_push_if_not_exists(void **array_ptr, const void *data)

---------
<br/>

#### array_reverse()


    void array_reverse(void *array)

---------
<br/>

#### key_value_get()


    char * key_value_get(const char *string, const char *key)

---------
<br/>

#### split()


    char ** split(char *string, char *delim)

---------
<br/>

#### getlinecount()


    long getlinecount(FILE *f)

---------
<br/>

#### backtrace()


    void backtrace()

---------
<br/>

<br/><br/><br/>

