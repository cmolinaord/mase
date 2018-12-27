# Isentropic compressible fluid dynamics simulation
# Computational Engineering 2018 - Computational Fluid Dynamics
# MASE - Master's degree in Space and Aeronautical Engineering
# Carlos Molina 2018

# Definition of the Domain
#      ________________________
# ny-1| ->                     |
#    ·| ->    N                |
#    ·| ->  W P E              |
#    3| ->    S                | H
#    2| ->                     |
#    1| ->                     |
#    0| ->_____________________|
#      0 1 2 3     · · ·   nx-1
#                 L

import time
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import Circle
import tools
import cfd
import const as c

# Input parameters
ITERMAX = 200
nx = 60
ny = 30
precission = 1e-3

args = sys.argv
if len(args) > 1:
	ITERMAX = int(args[1])
if len(args) > 3:
	nx = int(args[2])
	ny = int(args[3])
if len(args) > 4:
	precission = float(args[4])
elif len(args) > 5:
	print("ERROR: Too many input arguments")
	print("Exiting...")
	exit()

compressible = True

# Domain parameters (m)
L = 20
H = 10
# Solid objects
center = [L/2, H/2]
radius = min(L,H)/6
o = tools.obstacle(center, radius)

# Printing input data
print("ITERMAX = %i" % ITERMAX)
print("L = %1.1f" % L)
print("H = %1.1f" % H)
print("[nx,ny] = [%i,%i]" % (nx,ny))
print("Precission = %1.1E" % precission)
print("Compressible? %s" % ('Yes' if compressible else 'No'))
print("###########################")
print(" ")


tic = time.perf_counter()

# Inlet face
Vin = 3

# Creation of the World
w = tools.world(L, H, nx, ny)
# Creation of the fluid
f = tools.fluid(w, c.p0, c.T0, Vin)

# Boundary conditions computation
w, f = cfd.boundary(w, f, o)

# Iteration computing
iter = 0
error = 1

while error > precission and iter < ITERMAX:
	iter += 1
	for i in range(1, ny-1): # rows
		for j in range (1, nx-1): # columns
			f.phi_1[i,j] = tools.calc_phi(f.phi, f.rho, c.rho0, i, j, w)
	# Last column boundary condition (normal outflow)
	f.phi_1[:,nx-1] = f.phi[:,nx-2]
	#print(np.round(phi*100)/100)
	# Meassure error
	error = np.max(np.abs(f.phi_1 - f.phi))

	for i in range(1, ny-1): # rows
		for j in range (1, nx-1): # columns
			f.Vx[i,j], f.Vy[i,j] = tools.calc_vel(f.phi_1, f.rho, c.rho0, i, j, w, Vin)
			f.V[i,j] = np.sqrt(f.Vx[i,j]**2 + f.Vx[i,j]**2)
			# Energy conservation (calculated temperature)
			if compressible:
				f.T[i,j] = c.T0 + 0.5*(Vin**2 - f.V[i,j]**2)/c.c_p
				# Isentropic condition (pressure calculated)
				f.p[i,j] = c.p0 * (f.T[i,j]/c.T0)**c.gamma_exp
	# Compute new pressure (if compressible case)
	if compressible:
		f.rho = tools.density(f.p, f.T, c.R)
	f.phi[:,:] = f.phi_1[:,:]
	if iter % 10 == 0:
		print("Iteration %i: maximum error: %2.4e" %(iter, error))

toc = time.perf_counter() - tic
print("Elapsed time: %1.2fs" % toc)

cmap = plt.get_cmap('jet')

fig1 = plt.figure()
ax1 = fig1.gca()
ax1.streamplot(w.xv, w.yv, f.Vx, f.Vy, density=[0.5, 1])
im = ax1.pcolormesh(w.xv, w.yv, f.V, cmap=cmap)
fig1.colorbar(im, ax=ax1)
circ = Circle(center, radius, fill=False)
ax1.add_patch(circ)
ax1.axis("equal")
plt.title("Velocity (m/s)")

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.streamplot(w.xv, w.yv, f.Vx, f.Vy, density=[0.2, 0.5])
im = ax2.pcolormesh(w.xv, w.yv, f.rho, cmap=cmap)
fig2.colorbar(im, ax=ax2)
circ = Circle(center, radius, fill=False)
ax2.add_patch(circ)
ax2.axis("equal")
plt.title("Phi")

fig3 = plt.figure()
ax3 = fig3.gca()
ax3.plot(w.x,f.Vx[round(0.1*ny),:],'r-o')
ax3.plot(w.x,f.Vx[round(0.25*ny),:],'g-o')
ax3.plot(w.x,f.Vx[round(0.5*ny),:],'b-o')

plt.show()
