from __future__ import print_function
from __future__ import division

import sys
import time
from math import *
import numpy as np
import matplotlib.pyplot as plt

# import SIM5 library
sys.path.append('../../..')
import sim5

#parameters of the system
bh_mass = 10                    # BH mass [M_sun]
bh_spin = 0.99                  # BH spin
bh_dist = 1e3                   # BH distance [kpc]
bh_incl = 60.0                  # viewing angle (0=face-on, 90=edge-on)

disk_mdot = 0.1                 # accretion rate [Mdot_Eddington]
disk_alpha = 0.1                # alpha-viscosity

en_min = 1e-1                   # lower energy limit
en_max = 2e+1                   # upper energy limit
en_bins = 200                   # number of energy bins

# declare photon enerhy array in log scale
energies = np.logspace(log10(en_min), log10(en_max), en_bins)

# select a disk model (Novikove-Thorne in this case)
disk_model = sim5.DiskModel_ThinDisk(bh_mass, bh_spin, disk_mdot, disk_alpha)

# select emission model (blackbody in this case)
emission_model = sim5.DiskSpectrum_BlackBody()

# setup the raytracer for given disk and emission model
raytracer = sim5.DiskRaytrace(bh_mass, bh_spin, bh_dist, disk_model, emission_model)

# compute the spectrum
Ivf, Iv0 = raytracer.spectrum(bh_incl, energies, limbdk=1, flat=0)


# plot the result
plt.loglog(energies, Ivf, label='phenomenological hardening')
plt.loglog(energies, Iv0, label='without hardening')
plt.legend(frameon=False, loc='lower left')
plt.show()

