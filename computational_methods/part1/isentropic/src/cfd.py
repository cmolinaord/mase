import numpy as np
import const as c
import tools

def boundary(w, f, o):
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
			d = np.sqrt((w.x[j] - o.c[0])**2 + (w.y[i] - o.c[1])**2)
			if d < o.r:
				w.solid[i,j] = True

	return w, f

def gauss_seidel(w, f, opt):
	# Iteration computing
	iter = 0
	error = 1

	while error > opt.precission and iter < opt.itermax:
		iter += 1
		for i in range(1, w.ny-1): # rows
			for j in range (1, w.nx-1): # columns
				f.phi_1[i,j] = tools.calc_phi(f.phi, f.rho, c.rho0, i, j, w)
		# Last column boundary condition (normal outflow)
		f.phi_1[:,w.nx-1] = f.phi[:,w.nx-2]

		# Meassure error
		error = np.max(np.abs(f.phi_1 - f.phi))

		for i in range(1, w.ny-1): # rows
			for j in range (1, w.nx-1): # columns
				f.Vx[i,j], f.Vy[i,j] = tools.calc_vel(f.phi_1, f.rho, c.rho0, i, j, w, f.Vin)
				f.V[i,j] = np.sqrt(f.Vx[i,j]**2 + f.Vx[i,j]**2)

				if opt.compressible:
					# From energy conservation (calculated temperature)
					f.T[i,j] = c.T0 + 0.5*(f.Vin**2 - f.V[i,j]**2)/c.c_p
					# Isentropic condition (pressure calculated)
					f.p[i,j] = c.p0 * (f.T[i,j]/c.T0)**c.gamma_exp

		# Compute new pressure (if compressible)
		if opt.compressible:
			f.rho = tools.density(f.p, f.T, c.R)

		f.phi[:,:] = f.phi_1[:,:]
		if iter % 10 == 0 and opt.verbose == True:
			print("Iteration %i: maximum error: %2.4e" %(iter, error))
	return f
