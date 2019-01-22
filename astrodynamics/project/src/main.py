import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth, Jupiter, Sun
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter, plot
from astropy import time
from poliastro.util import norm
from common import *
from tools import Earth_to_Jupiter

#plt.ion()  # To immediately show plots
plt.style.use("seaborn")  # Recommended

now = time.Time.now()

from poliastro.maneuver import Maneuver

# Minimum Delta-V to reach Jupiter with Hohmann transfer
hoh = Maneuver.hohmann(ss_Earth, R_J)
dv_min = hoh[0][1]

# DeltaV to obtain Solar escape velocity
# (considering tangent burn)
escape_v = np.sqrt(2*Sun.k/R_E)
dv_escape = escape_v - norm(v_Earth)

# DeltaV of New Horizons mission
dv_NH_1 = 16.26*u.km/u.s

ss_trans, t, T, t_trans = Earth_to_Jupiter(dv_NH_1)

# Plotting
fig, ax = plt.subplots()
ax.grid(True)
ax.set_title("Transfer orbit from Earth to Jupiter with deltaV = %2.2f Km/s" % dv_NH_1.value)

op = OrbitPlotter(ax)
op.plot(ss_Earth, label="Now")
op.plot(ss_Jupiter, label="Now")
op.plot(ss_trans, label="Now")

op.plot(ss_Earth.propagate(t), label="At launch")
op.plot(ss_Jupiter.propagate(t), label="At launch")

op.plot(ss_Earth.propagate(t + t_trans), label="At arrival")
op.plot(ss_Jupiter.propagate(t + t_trans), label="At arrival")
op.plot(ss_trans.propagate(t_trans), label="At arrival")


plt.show()
