import numpy as np
from numpy import arccos
from astropy import units as u
from poliastro.util import norm
from common import *
from poliastro.maneuver import Maneuver
from poliastro.twobody.angles import nu_to_M

def angle(a,b):
	a = (a / norm(a))
	b = (b / norm(b))
	return arccos(np.dot(a,b))

def Earth_to_Jupiter(delta_v):
	"""
	Calculates the manuever and the launch windows from Earth to Jupiter,
	assuming tangential impulsive maneuver from Earth orbit around the Sun

	Earth_to_Jupiter(dv):

	    'dv' is given with velocity units (u.m/u.s)

	    Returns:
	        ss_transfer = trasfer orbit object
		  t = time to the first next launch window
		  T = period between possible launch windows
		  t_transfer = time of transfer from launch to Jupiter arrival
	"""
	dv = delta_v * v_Earth_dir

	# Initial transfer orbit calculation (to obtain transfer time)
	man_tmp = Maneuver.impulse(dv)
	transfer = ss_Earth.apply_maneuver(man_tmp)

	# Angles calculation
	# nu_Jupiter: nu of transfer orbit where it crosses Jupiter
	nu_Jupiter = (np.math.acos((transfer.p/R_J - 1)/transfer.ecc)*u.rad)
	ang_vel_Earth = (norm(v_Earth).value / R_E.value)*u.rad/u.s
	ang_vel_Jupiter = (norm(v_Jupiter).value / R_J.value)*u.rad/u.s

	# Transfer time
	M = nu_to_M(nu_Jupiter, transfer.ecc)
	t_transfer = M/transfer.n
	# Delta_nu at launch (Final nu - nu that Jupiter walks during transfer)
	delta_nu_at_launch = nu_Jupiter - ang_vel_Jupiter * t_transfer

	# Angles in radians
	nu_init = angle(r_Earth, r_Jupiter)*u.rad
	# Next launch window in t:
	t = (nu_init - delta_nu_at_launch) / (ang_vel_Earth - ang_vel_Jupiter)
	# Time between windows
	T = (2*np.pi*u.rad) / (ang_vel_Earth - ang_vel_Jupiter)

	# If time to the windows is negative (in the past), we should wait until next window
	while t < 0:
		t = t + T

	# Earth at launch
	ss_Earth_launch = ss_Earth.propagate(t)
	launch_dir = ss_Earth_launch.v / norm(ss_Earth_launch.v)
	dv = delta_v * launch_dir

	# Maneuver calculation
	man = Maneuver((t, dv))
	ss_transfer = ss_Earth.apply_maneuver(man)

	return ss_transfer, t, T, t_transfer
