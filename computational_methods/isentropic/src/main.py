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


import sys
import numpy as np
import tools
import matplotlib.pyplot as plt
from matplotlib.pyplot import Circle

# Input parameters
ITERMAX = 200
nx = 60
ny = 30
precission = 1e-3

args = sys.argv
if len(args) < 3:
	ITERMAX = int(args[1])
elif len(args) < 5:
	nx = int(args[2])
	ny = int(args[3])
elif len(args) < 6:
	precission = int(args[4])
else:
	print("ERROR: Too many input arguments")
	print("Exiting...")
	exit()

compressible = True

# Domain parameters (m)
L = 20
H = 10

# Physical constants
R = 287.058
c_p = 1.0035*1000 # Isobaric mass heat capacity (c_p) (J*kg^-1*K^-1)
gamma = 1.4 # Heat capacity ratio for dry air at room temperature
gamma_exp = gamma / (gamma-1)

# Physical input values
T0 = 273 # Initial temperature (Kelvin)
p0 = 101325 # Initial pressure (Pa)
rho0 = tools.density(p0, T0, R) # Intial density (kg/m^3)

# Inlet face
Vin = 3

# Creation of the World
w = tools.world(L, H, nx, ny)
phi 	= w.create_matrix(0) # Create a null matrix for phi
phi_new = w.create_matrix(0)
p 		= w.create_matrix(p0)
T 		= w.create_matrix(T0)
rho 	= w.create_matrix(rho0)
Vx 		= w.create_matrix(0)
Vy 		= w.create_matrix(0)
V 		= w.create_matrix(0)

# Initialization of boundary condition for Phi
# Inlet (West border)
for i in range(1, ny):
	phi[i,:] = phi[i-1,0] + Vin*w.dx
# North and South borders (Same phi as in the first column)
phi[0,1:] 	 = phi[0,0]
phi[ny-1,1:] = phi[ny-1,0]
# Copy initial conditions to phi_new
phi_new[:,:] = phi[:,:]

# Solid objects
center = [L/2, H/2]
radius = min(L,H)/6

w.solid[:,:] = False
# Defined walls in North and South
w.solid[0,:] = True
w.solid[ny-1,:] = True
# Defined cylinder in the middle
#Solid[5:45,20] = True
for i in range(1,ny-1):
	for j in range(1,nx-1):
		dist = np.sqrt((w.x[j]-center[0])**2 + (w.y[i]-center[1])**2)
		if dist < radius:
			w.solid[i,j] = True

# Iteration computing
iter = 0
error = 1

#fig2 = plt.figure()
#ax2 = fig2.gca()
#circ = Circle(center, radius, fill=False)
#ax2.add_patch(circ)
#ax2.axis("equal")
#plt.title("Density")


while error > precission and iter < ITERMAX:
	iter += 1
	for i in range(1, ny-1): # rows
		for j in range (1, nx-1): # columns
			phi_new[i,j] = tools.calc_phi(phi, rho, rho0, i, j, w)
	# Last column boundary condition (normal outflow)
	phi_new[:,nx-1] = phi[:,nx-2]
	#print(np.round(phi*100)/100)
	# Meassure error
	error = np.max(np.abs(phi_new - phi))

	for i in range(1, ny-1): # rows
		for j in range (1, nx-1): # columns
			Vx[i,j], Vy[i,j] = tools.calc_vel(phi_new, rho, rho0, i, j, w, Vin)
			V[i,j]	= np.sqrt(Vx[i,j]**2 + Vx[i,j]**2)
			# Energy conservation (calculated temperature)
			if compressible:
				T[i,j] = T0 + 0.5*(Vin**2 - V[i,j]**2)/c_p
				# Isentropic condition (pressure calculated)
				p[i,j] = p0 * (T[i,j]/T0)**gamma_exp
	# Compute new pressure (if compressible case)
	if compressible:
		rho = tools.density(p, T, R)
	phi[:,:] = phi_new[:,:]
	if iter % 10 == 0:
		print("Iteration %i: maximum error: %2.4e" %(iter, error))

	# Plots
	#ax2.streamplot(xv, yv, Vx, Vy, density=[0.5, 1])
	#im = ax2.pcolormesh(xv, yv, phi, cmap=cmap)
	#plt.pause(0.05)

cmap = plt.get_cmap('PiYG')

fig1 = plt.figure()
ax1 = fig1.gca()
ax1.streamplot(w.xv, w.yv, Vx, Vy, density=[0.5, 1])
im = ax1.pcolormesh(w.xv, w.yv, V, cmap=cmap)
fig1.colorbar(im, ax=ax1)
circ = Circle(center, radius, fill=False)
ax1.add_patch(circ)
ax1.axis("equal")
plt.title("Velocity (m/s)")

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.streamplot(w.xv, w.yv, Vx, Vy, density=[0.5, 1])
im = ax2.pcolormesh(w.xv, w.yv, phi, cmap=cmap)
fig2.colorbar(im, ax=ax2)
circ = Circle(center, radius, fill=False)
ax2.add_patch(circ)
ax2.axis("equal")
plt.title("Phi")


fig3 = plt.figure()
ax3 = fig3.gca()
ax3.plot(w.x,Vx[round(0.1*ny),:],'r-o')
ax3.plot(w.x,Vx[round(0.25*ny),:],'g-o')
ax3.plot(w.x,Vx[round(0.5*ny),:],'b-o')

plt.show()
