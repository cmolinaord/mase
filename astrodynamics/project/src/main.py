# Computation of possible future launch windows
# to reach Jupiter from Earth using Poliastro
#   Carlos Molina (2019)
#   carlosmolina.ord@gmail.com
#   https://github.com/cmolinaord/mase/tree/master/astrodynamics/project/src
#
# Universitat Politecnica de Catalunya
# Master's Degree in Space and Aeronautical Engineering (MASE)
# Astrodynamics final project
#
# Packages needed to run:
# numpy
# matplotlib
# astropy
# poliastro

import numpy as np
import matplotlib.pyplot as plt
from astropy import time
from astropy import units as u
from poliastro.util import norm
from poliastro.bodies import Earth, Jupiter, Sun
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter, plot
from poliastro.maneuver import Maneuver
from common import *
from tools import Earth_to_Jupiter

plt.style.use("seaborn")  # Recommended

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

print("Time from now to the next launch window = %2.3f d" % t.to(u.d).value)
print("Time between launch windows = %2.3f d" % T.to(u.d).value)
print("Transfer time from Earth to Jupiter = %2.3f d" % t_trans.to(u.d).value)

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

# Second plot, for different delta_V values

dv = np.linspace(norm(dv_min), 18*u.km/u.s, num=9)
tt = np.zeros([len(dv),1]) # Time to next launch window
TT = np.zeros([len(dv),1]) # Transfer time
total_time = np.zeros([1,len(dv)]) # Total time needed

for i in range(len(dv)):
	ss_trans, tt[i], T, TT[i] = Earth_to_Jupiter(dv[i])

tt = (tt*u.s).to(u.d).value
TT = (TT*u.s).to(u.d).value

total_time = tt + TT

fig2, ax2 = plt.subplots()
ax2.grid(True)
ax2.set_title("Transfers times as function of deltaV")

ax2.plot(dv, tt, "-o", label="Time to next launch window")
ax2.plot(dv, TT, "-o", label="Transfer time")
ax2.plot(dv, total_time, "-o", label="Total time from today")
ax2.legend()
ax2.grid(True)
ax2.set_xlabel("Delta_V (km/s)")
ax2.set_ylabel("Time (d)")
plt.show()
