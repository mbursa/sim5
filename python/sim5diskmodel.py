# Basic Disk models.
#
# Module defines few basic classes for accretion disk models:
#  - an empty class for further inheritance
#  - class for relativistic thin disk model (Novikov-Thorne)
#  - class for an external model loaded from a library
# 
# This file is a part of SIM5 library. 
# See README and LICENCE file for details.

import sys
import ctypes
import importlib
import logging
import sim5lib as sim5



class DiskModel:
    """
    The base class for deriving disk models.

    DiskModel class defined all basic methods needed for a disk model such as methods
    for giving flux, surface density, angular momentum, radial velocity, vertical 
    scale-height, and its radial derivative. Those methods are used by other classes 
    (e.g. in raytracing procedure) to extract disk properties at certain radii.
    
    In all methods that take radius (R) as an argument, this radius shall by convention
    be interpreted as 'equatorial radius', i.e. the radial Boyer-Lindquist coordinate
    measured in the equatorial plane (at theta=pi/2).
    
    All methods/attributes giving physical quantities shall use CGS unit system or 
    return a dimensionless number.
    
    Attributes:
        name: Name of the model (string).
    """

    def __init__(self):
        self.name = '(disk model base class)'
    #end of def


    # radiative flux [erg/s/cm2]
    def flux(self, R): return 0.0

    # effective temperature [K]
    def t_eff(self, R): return (self.flux(R)/5.670400e-05)**0.25  # T_eff = (F/sigma_sb)^1/4

    # surface density [g/cm2]
    def sigma(self, R): return 0.0

    # specific angular momentum [???]
    def l(self, R): return 0.0

    # radial velocity [c]
    def vr(self, R): return 0.0

    # disk scale-height [rg]
    def h(self, R): return 0.0

    # derivative (slope) of disk surface
    def dhdr(self, R): return 0.0

#end class




class DiskModel_ThinDisk(DiskModel):
    """
    Disk model that links an external library/module.
    """


    def __init__(self, bh_mass, bh_spin, mdot_L, alpha, options=0):
        sim5.disk_nt_setup(bh_mass, bh_spin, mdot_L, alpha, options or 0)
        self.name  = 'Novikov-Thorne'
        self.mdot  = sim5.disk_nt_mdot()
        self.lumi  = sim5.disk_nt_lumi()
        self.r_min = sim5.disk_nt_r_min()
        logging.info("Disk model: '%s' M=%.1f a=%.3f mdot=%.3e (%.5e g/s) lum=%.3e", self.name, bh_mass, bh_spin, self.mdot, self.mdot*bh_mass*sim5.Mdot_Edd, self.lumi)
    #end of def

    def flux(self, R): return sim5.disk_nt_flux(R)

    def sigma(self, R): return sim5.disk_nt_sigma(R)

    def l(self, R): return sim5.disk_nt_ell(R)

    def vr(self, R): return 0.0

    def h(self, R): return 0.0

    def dhdr(self, R): return 0.0
#end class




class DiskModel_External:
    """
    Disk model that links an external library/module.
    """

    def __init__(self, model_lib, bh_mass, bh_spin, options=''):
        """
        Initializes an external disk model.

        Args:
            model_lib: a shared library or a python module that provides disk model methods
            bh_mass: black hole mass parameter [M_sun]
            bh_spin: black hole spin parameter [0..1]
            options: other options tht are passed to the model library (a sequence of key=value pairs delimited by comma)
        """
        self.lib = None

        try:
            if (model_lib and model_lib.endswith('.so')):
                # load shared library and setup c_types interface for its functions
                # example: http://gestaltrevision.be/wiki/python/cfun
                self.lib = ctypes.cdll.LoadLibrary(model_lib)

                self.lib.diskmodel_init.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_char_p]
                self.lib.diskmodel_init.restype = None

                self.lib.diskmodel_done.argtypes = []
                self.lib.diskmodel_done.restype = None

                self.lib.diskmodel_name.argtypes = []
                self.lib.diskmodel_name.restype = ctypes.c_char_p

                self.lib.diskmodel_r_min.argtypes = []
                self.lib.diskmodel_r_min.restype = ctypes.c_double

                self.lib.diskmodel_mdot.argtypes = []
                self.lib.diskmodel_mdot.restype = ctypes.c_double

                self.lib.diskmodel_lumi.argtypes = []
                self.lib.diskmodel_lumi.restype = ctypes.c_double

                self.lib.diskmodel_flux.argtypes = [ctypes.c_double]
                self.lib.diskmodel_flux.restype = ctypes.c_double

                self.lib.diskmodel_sigma.argtypes = [ctypes.c_double]
                self.lib.diskmodel_sigma.restype = ctypes.c_double

                self.lib.diskmodel_l.argtypes = [ctypes.c_double]
                self.lib.diskmodel_l.restype = ctypes.c_double

                self.lib.diskmodel_vr.argtypes = [ctypes.c_double]
                self.lib.diskmodel_vr.restype = ctypes.c_double

                self.lib.diskmodel_h.argtypes = [ctypes.c_double]
                self.lib.diskmodel_h.restype = ctypes.c_double

                self.lib.diskmodel_dhdr.argtypes = [ctypes.c_double]
                self.lib.diskmodel_dhdr.restype = ctypes.c_double
            elif (model_lib and model_lib.endswith('.py')):
                self.lib = importlib.import_module(model_lib)
            #end if                
        except AttributeError as e:
            print('Missing required function in disk module.\n', e.args)
            raise

        # init model
        if (self.lib):
            self.lib.diskmodel_init(bh_mass, bh_spin, options)
            self.r_min = self.lib.diskmodel_r_min()
            self.name  = self.lib.diskmodel_name()
            self.mdot  = self.lib.diskmodel_mdot()
            self.lumi  = self.lib.diskmodel_lumi()
            logging.info("# Disk model: '%s' M=%.1f a=%.3f mdot=%.3e (%.5e) lum=%.3e", self.name, bh_mass, bh_spin, self.mdot, self.mdot*bh_mass*2.225475942e+18, self.lumi)
        #end if

    #end of def


    def __del__(self):
        if (self.lib and self.lib.diskmodel_done): self.lib.diskmodel_done()
    #end of def

    def close(self):
        if (self.lib and self.lib.diskmodel_done): self.lib.diskmodel_done()
        self.lib = None
    #end of def

    def flux(self, R): return self.lib.diskmodel_flux(R) if self.lib else 0.0

    def sigma(self, R): return self.lib.diskmodel_sigma(R) if self.lib else 0.0

    def l(self, R): return self.lib.diskmodel_l(R) if self.lib else 0.0

    def vr(self, R): return self.lib.diskmodel_vr(R) if self.lib else 0.0

    def h(self, R): return self.lib.diskmodel_h(R) if self.lib else 0.0

    def dhdr(self, R): return self.lib.diskmodel_dhdr(R) if self.lib else 0.0

#end of class





'''
# function to test the class
def test():
    if (len(sys.argv) < 4):
        print "Usage: %s sd|nt spin mdot"
        sys.exit()
    #end if

    model = sys.argv[1]    
    spin  = float(sys.argv[2])
    mdot  = float(sys.argv[3])
    alpha = 0.1

    m = DiskModel('./flux-'+model+'.so', 10, 0.0, 'alpha=%f'%(alpha))

    print '# model dump:', m.name
    print '# rmin:', m.r_min
    print '# mdot:', m.mdot
    print '# lumi:', m.lumi
    print '# alpha:', alpha
    print '# [radius  flux  sigma  ell  vr  h  dhdr]'
    r = m.r_min
    while (r < 1e4):
        print r, m.flux(r), m.sigma(r), m.l(r), m.vr(r), m.h(r), m.dhdr(r)
        r *= 1.05
    #end while
    

    print 
#end of def


# if executed as main file then
# run the main function and exit with returned error code
if __name__ == "__main__": sys.exit(__test())
'''


