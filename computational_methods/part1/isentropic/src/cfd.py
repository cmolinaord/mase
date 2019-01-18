import numpy as np
import const as c
from tools import *
from copy import deepcopy

def boundary(w, f, obs):
	# Initialization of boundary condition for Phi
	# Inlet (West border)
	for i in range(1, w.ny):
		f.phi[i,:] 	= f.phi[i-1,0] + f.Vin*w.dx
	# North and South borders (Same phi as in the first column)
	f.phi[0,1:] 	= f.phi[0,0]
	f.phi[w.ny-1,1:] 	= f.phi[w.ny-1,0]
	# Copy initial conditions to phi_1
	f.phi_1[:,:] 	= f.phi[:,:]

	# Solid computation
	w.solid[:,:] = False
	# Defined walls in North and South
	w.solid[0,:] = True
	w.solid[w.ny-1,:] = True

	for i in range(1,w.ny-1):
		for j in range(1,w.nx-1):
			d = np.sqrt((w.x[j] - obs.c[0])**2 + (w.y[i] - obs.c[1])**2)
			if d < obs.r:
				w.solid[i,j] = True

	return w, f

def gauss_seidel(w, f, opt):
	# Iteration computing
	iter = 0
	error = 1

	while error > opt.precission and iter < opt.itermax:
		iter += 1
		# Phi calc
		for i in range(1, w.ny-1): # rows
			for j in range(1, w.nx-1): # columns
				f = calc_phi(f, w, c, i, j, w)
		# Last column boundary condition (normal outflow)
		f.phi_1[:,w.nx-1] = f.phi[:,w.nx-2]

		# Meassure error
		error = np.max(np.abs(f.phi_1 - f.phi))
		# Velocity calc
		for i in range(1, w.ny-1): # rows
			for j in range (1, w.nx-1): # columns
				f = calc_vel(f, c, i, j, w)

				if opt.compressible:
					# From energy conservation (calculated temperature)
					f.T[i,j] = c.T0 + 0.5*(f.Vin**2 - f.V[i,j]**2)/c.c_p
					# Isentropic condition (pressure calculated)
					f.p[i,j] = c.p0 * (f.T[i,j]/c.T0)**c.gamma_exp

		# Compute new pressure (if compressible)
		if opt.compressible:
			f.rho = density(f.p, f.T, c.R)

		f.phi[:,:] = f.phi_1[:,:]
		if iter % 10 == 0 and opt.verbose == True:
			print("  Iteration %i: maximum error: %2.4e" %(iter, error))
	res = results()
	res.iters = iter
	res.error = error
	return f, res

def solid_circulation(f, w):
	circ = 0
	c = 0
	# Sum up the circulation value
	# border velocities * lenght (clockwise sign criteria)
	c = f.vn*w.dx - f.ve*w.dy - f.vs*w.dx + f.vw*w.dy
	circ = np.sum(c*w.solid)
	return circ

def circulation_computation(obstacle, options):
	# This function computes the whole domain, with a reduced
	# resolution, and calculates the circulation around the obstacle
	# for different test phi_solid inside the obstacle

	N = 4 # Number of tries
	rX = 0.8 # Reduction of mesh resolution
	pX = 0.6 # Reduction in precission
	iX = 0.6 # Reduction in iterations

	# Create local copy of this objects, to not affect global ones
	obs = deepcopy(obstacle)
	opt = deepcopy(options)

	opt.nx 		= int(opt.nx * rX)
	opt.ny 		= int(opt.ny * rX)
	opt.precission 	= opt.precission / pX
	opt.itermax 	= int(opt.itermax * iX)
	if opt.nx*opt.ny*opt.itermax < 100000:
		opt.verbose 	= False

	print("Starting circulation computation with:")
	print("    %i tries of:" % N)
	print("    %1.1f reduced mesh resolution: (%i,%i)" % (rX,opt.nx,opt.ny))
	print("    %1.1f reduced precision: %1.2e" % (pX,opt.precission))
	print("    %1.1f reduced number of iterations: %i" % (iX,opt.itermax))
	print(" ")
	print("Circulation wanted: %1.2f" % obs.circ)

	# Phi values to try (from 0 to Vin*dx*ny)
	phi_try = np.linspace(0, opt.Vin*opt.L/opt.nx*opt.ny, N)
	circ = np.zeros([N])

	# Common world
	w_aux = world(opt)
	f = fluid(w_aux, c, opt)
	w_aux, f = boundary(w_aux, f, obs)

	for i in range(N):
		w_aux.phi_solid = phi_try[i]
		f, res = gauss_seidel(w_aux, f, opt)
		circ[i] = solid_circulation(f, w_aux)
		print("Try(%i) iters: %i" % (i,res.iters))

	m = (circ[-1] - circ[0])/(phi_try[-1] - phi_try[0])
	n = circ[0] - m*phi_try[0]

	phi_solid = (obs.circ - n) / m
	print("Computed phi_solid = %1.2f\n" % phi_solid)
	return phi_solid
