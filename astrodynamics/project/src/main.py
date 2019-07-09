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
from poliastro.bodies import Earth, Jupiter, Sun, Mars
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter, plot
from poliastro.maneuver import Maneuver
from tools import Tangent_burn

plt.style.use("seaborn")  # Recommended

A = Earth
B = Mars

ss_A = Orbit.from_body_ephem(A)
ss_B = Orbit.from_body_ephem(B)

# Minimum Delta-V to reach Jupiter with Hohmann transfer
hoh = Maneuver.hohmann(ss_A, norm(ss_B.r.to(u.m)))
dv_min = hoh[0][1]

# DeltaV of New Horizons mission
dv_NH_1 = 16.26*u.km/u.s
dv_NH_1 = 3.8*u.km/u.s

ss, t, T, t_trans = Tangent_burn(A, B, dv_NH_1)

print("Time from now to the next launch window = %2.3f d" % t.to(u.d).value)
print("Time between launch windows = %2.3f d" % T.to(u.d).value)
print("Transfer time from %s to %s = %2.3f d" % (A.name, B.name, t_trans.to(u.d).value))

# Plotting
fig, ax = plt.subplots()
ax.grid(True)
ax.set_title("Transfer orbit from %s to %s with deltaV = %2.2f Km/s" % (A.name, B.name, dv_NH_1.value))

op = OrbitPlotter(ax)
op.plot(ss.A, label="Now")
op.plot(ss.B, label="Now")

op.plot(ss.A.propagate(t), label="At launch")
op.plot(ss.B.propagate(t), label="At launch")
op.plot(ss.trans, label="At launch")

op.plot(ss.A.propagate(t + t_trans), label="At arrival")
op.plot(ss.B.propagate(t + t_trans), label="At arrival")
op.plot(ss.trans.propagate(t_trans), label="At arrival")

plt.show()

# Second plot, for different delta_V values

dv = np.linspace(norm(dv_min), 18*u.km/u.s, num=9)
tt = np.zeros([len(dv),1]) # Time to next launch window
TT = np.zeros([len(dv),1]) # Transfer time
total_time = np.zeros([1,len(dv)]) # Total time needed

for i in range(len(dv)):
	ss, tt[i], T, TT[i] = Tangent_burn(A, B, dv[i])

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
