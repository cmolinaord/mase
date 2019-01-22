import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from poliastro.bodies import Earth, Jupiter, Sun
from poliastro.twobody import Orbit
from poliastro.plotting import OrbitPlotter, plot
from poliastro import iod
from astropy import time
from poliastro.util import norm


plt.ion()  # To immediately show plots
plt.style.use("seaborn")  # Recommended


epoch = time.Time("2015-05-09 10:43")

now = time.Time.now()


ss_Earth = Orbit.from_body_ephem(Earth)
ss_Jupiter = Orbit.from_body_ephem(Jupiter)

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
op.plot(ss_exit)

#for v in [9, 10, 11, 12, 13, 14, 15]:
#	dv = v * v_tg * u.km / u.s
#	ss_Earth.apply_maneuver(Maneuver.impulse(dv))
#	op.plot()
