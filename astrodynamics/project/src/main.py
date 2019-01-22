import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth, Jupiter, Sun
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter, plot
from poliastro import iod
from astropy import time
from poliastro.util import norm
from math import acos


plt.ion()  # To immediately show plots
plt.style.use("seaborn")  # Recommended


epoch = time.Time("2015-05-09 10:43")

now = time.Time.now()

ss_Earth = Orbit.from_body_ephem(Earth)
r_Earth = ss_Earth.r.to(u.m)
R_E = norm(r_Earth)
v_Earth = ss_Earth.v.to(u.m/u.s)
v_Earth_dir = v_Earth / norm(v_Earth)

ss_Jupiter = Orbit.from_body_ephem(Jupiter)
r_Jupiter = ss_Jupiter.r.to(u.m)
R_J = norm(r_Jupiter)
v_Jupiter = ss_Jupiter.v.to(u.m/u.s)

# Ultima Thule orbit elements at epoch 26 April 2019 00:00 UTC

from UltimaThule import orbit as ss_MU69

op = OrbitPlotter()
op.plot(ss_Earth)
op.plot(ss_Jupiter)
# op.plot(ss_MU69)

print("Earth: ", ss_Earth.state.nu.to(u.deg))
print("Jupiter: ", ss_Jupiter.state.nu.to(u.deg))
print("Ultima Thule: ", ss_MU69.state.nu.to(u.deg))

from poliastro.maneuver import Maneuver

# Create Hohmann maneuver to calculate minimum deltaV to reach Jupiter
hoh = Maneuver.hohmann(ss_Earth, R_jupiter)
dv = hoh[0][1]
print("Minimum delta-v to reach Jupiter with Hohmann = " % norm(dv))

escape_v = np.sqrt(2*Sun.k/r_E)

# DeltaV escape nedded to actually escape Earth
# (considering tangent burn)
dv_escape = escape_v - norm(v_Earth)

dv_NH_1 = 16.26*u.km/u.s

dv = dv_NH_1 * v_Earth_dir
# First calculation of time to Jupiter using New Horizons escape velocity
man_NH_1 = Maneuver.impulse(dv)
ss_NH_1 = ss_Earth.apply_maneuver(man_NH_1)
# nu at Jupiter orbit cross
nu_at_Jupiter = (np.math.acos((ss_NH_1.state.p/R_jupiter - 1)/ss_NH_1.state.ecc)*u.rad)

ang_vel_Earth = (norm(v_Earth).value / R_E.value)*u.rad/u.s
ang_vel_Jupiter = (norm(v_Jupiter).value / R_J.value)*u.rad/u.s

# Angles in radians
nu_init = tools.angle(r_Earth, r_Jupiter)*u.rad
# Next launch window in t:
t = -(nu_at_Jupiter - nu_init) / (ang_vel_Earth - ang_vel_Jupiter)

T_window = (2*np.pi*u.rad) / (ang_vel_Earth - ang_vel_Jupiter)




op.plot(ss_exit)

#for v in [9, 10, 11, 12, 13, 14, 15]:
#	dv = v * v_tg * u.km / u.s
#	ss_Earth.apply_maneuver(Maneuver.impulse(dv))
#	op.plot()
