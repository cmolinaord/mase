import numpy as np
#from tools import object as o

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
