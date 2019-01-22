import numpy as np
from astropy import units as u
from poliastro.twobody import Orbit
from poliastro.twobody.angles import M_to_nu
from poliastro.bodies import Sun


MU69_ecc = 0.05046617*u.one	# Eccentricity
MU69_i = 2.44996247*u.deg	# Inclination
MU69_a = 6676495970*u.km 	# Semimajor axisº
MU69_raan = 159.04395317345*u.deg	# Right ascension of the ascending node
MU69_argp = 181.14813924539*u.deg	# Argument of perihelion
MU69_M = 310.91947579009*u.deg	# Mean anomaly
MU69_nu = M_to_nu(MU69_M, MU69_ecc)	# True anomaly

orbit = Orbit.from_classical(Sun, MU69_a, MU69_ecc, MU69_i, MU69_raan, MU69_argp, MU69_nu)
