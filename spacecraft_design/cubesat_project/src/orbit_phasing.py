import numpy as np
import matplotlib.pyplot as plt
plt.ion()  # To immediately show plots
from astropy import units as u
from poliastro.bodies import Earth, Mars, Sun
from poliastro.twobody import Orbit
plt.style.use("seaborn")  # Recommended
h = 625*u.km
r_circ = Earth.R + h
circ = Orbit.circular(Earth,h)
T_circ = circ.state.period
N = 9
delay = T_circ/N
T_ph = T_circ + delay
a_ph = np.power(T_ph/2/np.pi * np.sqrt(Earth.k),(2/3))
rp=Earth.R + h
ra = 2*a_ph - rp
e_ph = (ra-rp)/(ra+rp)
phasing = Orbit.from_classical(Earth,a_ph,e_ph,0*u.deg,0*u.deg,0*u.deg,0*u.deg)
from poliastro.plotting import OrbitPlotter
op = OrbitPlotter()
op.plot(circ, label="Final circular orbit")
op.plot(phasing, label="Phasing orbit")
