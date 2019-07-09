import numpy as np
from numpy import arccos
from astropy import units as u
from poliastro.util import norm
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver
from poliastro.twobody.angles import nu_to_M

class Orbits(object):
	def __init__(self, planetA, planetB):
		self.A = Orbit.from_body_ephem(planetA)
		self.B = Orbit.from_body_ephem(planetB)
		# Radius vector
		self.r_A = self.A.r.to(u.m)
		self.r_B = self.B.r.to(u.m)
		# Radius vector norm
		self.R_A = norm(self.r_A)
		self.R_B = norm(self.r_B)
		# Velocity vectors
		self.v_A = self.A.v.to(u.m/u.s)
		self.v_B = self.B.v.to(u.m/u.s)
		# Velocity vectors direction
		self.v_A_dir = self.v_A / norm(self.v_A)

def angle(a,b):
	a = (a / norm(a))
	b = (b / norm(b))
	return arccos(np.dot(a,b))

def Tangent_burn(planetA, planetB, delta_v):
	"""
	Calculates the manuever and the launch windows from planet A to planet B,
	assuming tangential impulsive maneuver from A orbit around the Sun

	tangent_burn(planetA, plnaet_B, dv):

	    'planet_X': are the planets origin and target
	    'dv' is given with velocity units (u.m/u.s)

	    Returns:
	        ss_transfer = trasfer orbit object
		  t = time to the first next launch window
		  T = period between possible launch windows
		  t_transfer = time of transfer from launch to arrival at B
	"""

	ss = Orbits(planetA, planetB)

	dv = delta_v * ss.v_A_dir

	# Initial transfer orbit calculation (to obtain transfer time)
	man_tmp = Maneuver.impulse(dv)
	transfer = ss.A.apply_maneuver(man_tmp)

	# Angles calculation
	# nu_B: nu of transfer orbit where it crosses B
	nu_B = (np.math.acos((transfer.p/ss.R_B - 1)/transfer.ecc)*u.rad)
	ang_vel_A = (norm(ss.v_A).value / ss.R_A.value)*u.rad/u.s
	ang_vel_B = (norm(ss.v_B).value / ss.R_B.value)*u.rad/u.s

	# Transfer time
	M = nu_to_M(nu_B, transfer.ecc)
	t_transfer = M/transfer.n
	# Delta_nu at launch (Final nu - nu that B walks during transfer)
	delta_nu_at_launch = nu_B - ang_vel_B * t_transfer

	# Angles in radians
	nu_init = angle(ss.r_A, ss.r_B)*u.rad
	# Next launch window in t:
	t = (nu_init - delta_nu_at_launch) / (ang_vel_A - ang_vel_B)
	# Time between windows
	T = (2*np.pi*u.rad) / (ang_vel_A - ang_vel_B)

	# If time to the windows is negative (in the past), we should wait until next window
	while t < 0:
		t = t + T

	# A at launch
	ss_A_launch = ss.A.propagate(t)
	launch_dir = ss_A_launch.v / norm(ss_A_launch.v)
	dv_launch = delta_v * launch_dir

	# Maneuver calculation
	man = Maneuver.impulse(dv_launch)
	ss_transfer = ss_A_launch.apply_maneuver(man)

	ss.trans = ss_transfer

	return ss, t, T, t_transfer
