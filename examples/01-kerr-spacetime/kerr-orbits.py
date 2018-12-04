from __future__ import division

import sys
from math import *
import numpy as np
import matplotlib.pyplot as plt

# import SIM5 library
sys.path.append('../../..')
import sim5

# make array of spin values from 0 to 1
bh_spin = np.linspace(0, 1, 200)

# compute 
# - black-hole event horizon radius,
# - marginally stable orbit radius,
# - marginally bound radius,
# - photon orbit radius
# for all spins in bh_spin array
rbh = map(sim5.r_bh, bh_spin)
rms = map(sim5.r_ms, bh_spin)
rmb = map(sim5.r_mb, bh_spin)
rph = map(sim5.r_ph, bh_spin)

# make a plot of rbh and rms as a function of spin
plt.title("Kerr orbits")
plt.plot(bh_spin, rms, label='marg. stable orbit')
plt.plot(bh_spin, rmb, label='marg. bound orbit')
plt.plot(bh_spin, rph, label='photon orbit')
plt.plot(bh_spin, rbh, label='event horizon')
plt.xlabel(r'BH spin')
plt.ylabel(r'radius [$r_{\rm g}$]')
plt.legend()

plt.show()

