# Spectral models for local emission of accretion disks.
#
# Module provides classes that define spectral model for radiation
# emerging from an accretion disk.
# 
# This file is a part of SIM5 library. 
# See README and LICENCE file for details.


from __future__ import division
import sys
import numpy as np



class DiskSpectrum:
    def __init__(self):
        pass
    #end of def


    def spectrum(self, T, m, f, E):
        """
        Computes spectrum for specified parameters by interpolating TLUSTY
        spectra on given energy grid.

        Makes interpolation on the parameter grid and returns weighted spectrum
        projected onto a given energy grid `E` for emission angle `m`.

        If [T,S,Q] point lies outside the spectral grid, a failover black-body spectrum
        is returned with hardening factor f.

        Args:
            T: temperature [K]
            m: cosine of emission angle
            f: hardening factor to be used in case of black-body failover
            E: array of energies to which to project the spectrum [keV]
        Returns:
            specific intensity grid [erg/cm2/s/keV/srad]
        """

        return np.zeros(len(E))
    #end of def
#end of class



class DiskSpectrum_BlackBody(DiskSpectrum):
    def __init__(self):
        pass
    #end of def


    def spectrum(self, T, m, f, E):
        """
        Computes spectrum for specified parameters by evaluating Planck formula
        for given temperature.

        Returns specific intensity Iv of a black body using Planck formula. The spectrum
        is modifies by a hardening factor, which redistributes photons from softer
        to higher energies while keeping the total flux, and by limb-darkening effect,
        which redistributed photons between emission angles (more emission along surface normal).

        Args:
            T: temperature [K]
            m: cosine of emission angle (if m<0 then isotropic emission is assumed)
            f: hardening factor
            E: array of energies to which to project the spectrum [keV]
        Returns:
            specific intensity [erg/cm2/s/keV/srad]
        """

        planck_h        = 6.626069e-27               # Planck constant [erg.s]
        kev2freq        = 2.417990e+17               # keV->Hz (1e3*electronvolt/planck_h)
        speed_of_light2 = 8.987554e+20               # square of speed of light [cm^2/s^2]
        boltzmann_k     = 1.380650e-16               # Boltzmann constant [erg/K]

        if (T < 1e2): return np.zeros(len(E))
        # make spectrum using black-body formula with hardening factor & limb darkening
        limbf = 0.5+0.75*m if (m>=0.0) else 1.0

        # calc Planck spectrum in units [erg/cm^2/s/keV/srad]
        with np.errstate(over='ignore'):
            Iv = limbf*2.0*planck_h*(kev2freq*E)**3/speed_of_light2/(f**4) * 1./(np.exp((planck_h*kev2freq*E)/(boltzmann_k*f*T))-1.0) * kev2freq

        return Iv
    #end of def
#end of class




