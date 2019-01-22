import numpy as np
from astropy import units as u
from poliastro.bodies import Earth, Jupiter, Sun
from poliastro.twobody import Orbit
from poliastro.util import norm

# Orbits
ss_Earth = Orbit.from_body_ephem(Earth)
ss_Jupiter = Orbit.from_body_ephem(Jupiter)

# Radius vector
r_Earth = ss_Earth.r.to(u.m)
r_Jupiter = ss_Jupiter.r.to(u.m)

# Radius vector norm
R_E = norm(r_Earth)
R_J = norm(r_Jupiter)

# Velocity vectors
v_Earth = ss_Earth.v.to(u.m/u.s)
v_Jupiter = ss_Jupiter.v.to(u.m/u.s)

# Velocity vectors direction
v_Earth_dir = v_Earth / norm(v_Earth)
