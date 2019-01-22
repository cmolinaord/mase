import numpy as np
from math import acos
from astropy import units as u
from poliastro.util import norm


def angle(a,b):
	a = (a / norm(a))
	b = (b / norm(b))
	return np.arccos(np.dot(a,b))
