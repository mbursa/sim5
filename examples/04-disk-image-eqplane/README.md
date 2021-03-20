
# Example 4 - Accretion disk image (equatorial plane)

This example illustrates how SIM5 library can be used in a C program and in a Python script to compute an image of an axially symmetric thin accretion disk around a black hole. The disk assumed to be located at the equatroial plane and to have zero geometrical thickness (an usual astrophysical assumption for disks accreting at low rates).
The equatorial plane setup allows a fast evaluation (up to ~10<sup>5</sup> photons per second) of the image because the radius from which the photons originate may be calculated by an inversion of an elliptic positional integral. Thus, by using 
```C
geodesic_init_inf(inclination, spin, alpha, beta, &gd, &error);
P = geodesic_find_midplane_crossing(&gd, 0);
r = geodesic_position_rad(&gd, P);
```
one can quickly and efficiently find the original radius of a photon coming out of the disk surface that has  impact parameters *alpha* and *beta*.

## Running in C
```bash
make
./disk-image 0.5 60 > image.dat
```

First argument is the black-hole spin (0 to 0.999), second argument is the angle of the line of sight of an observer with respect to symmetry axis (0 to 90 degrees).

The program computes the disk image and dumps the image data to standard output. Stdout can be redirected to a file. The format of each printed line is:  
`y  x  flux  g_factor`
where *x*, *y* are pixel coordinates, *flux* is a frequency-integrated radiative flux received a pixels and *g_factor* is the value of the energy shift due to relativistic effects.

## Running in Python
```bash
python disk-image.py 0.5 80
```

First argument is the black-hole spin (0 to 0.999), second argument is the line-of-sight angle (0 to 90 degrees). The script diaplays the plot of flux and g-factor.

## Extending the example
In addition to `geodesic_position_rad()` function, other positional parameters of the photon can be determined by `geodesic_position_azm()` (azimuthal angle) for non-axisymmteric applications, and `geodesic_timedelay()` for timing applications. Direction of the outgoing photon can be obtained from `geodesic_momentum()`.

