# Isentropic compressible fluid dynamics simulation
# Computational Engineering 2018 - Computational Fluid Dynamics
# MASE - Master's degree in Space and Aeronautical Engineering
# (C) Carlos Molina 2018

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

args = sys.argv
ITERMAX = int(args[1])

compressible = False

# Computational parameters
precission = 1e-3
nx = 120
ny = 60
iter = 0

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

# Creation of domain
phi, dx, dy = tools.create_mesh(L, H, nx, ny) # Create a null matrix for phi
phi 	= np.zeros([ny,nx])
phi_new = np.zeros([ny,nx])
p 	= np.zeros([ny,nx]) + p0
T 	= np.zeros([ny,nx]) + T0
rho 	= np.zeros([ny,nx]) + rho0
Vx 	= np.zeros([ny,nx])
Vy 	= np.zeros([ny,nx])
V 	= np.zeros([ny,nx])

# Initialization of boundary condition for Phi
# Inlet (West border)
for i in range(1, ny):
	phi[i,:] = phi[i-1,0] + Vin*dx
# North and South borders (Same phi as in the first column)
phi[0,1:] 	 = phi[0,0]
phi[ny-1,1:] = phi[ny-1,0]
# Copy initial conditions to phi_new
phi_new[:,:] = phi[:,:]

x = np.linspace(0, L, nx)
y = np.linspace(0, H, ny)
xv, yv = np.meshgrid(x, y)

# Solid objects
center = [L/2, H/2]
radius = min(L,H)/6

Solid, dx, dy = tools.create_mesh(L, H, nx, ny)
Solid[:,:] = False
# Defined walls in North and South
Solid[0,:] = True
Solid[ny-1,:] = True
# Defined cylinder in the middle
#Solid[5:45,20] = True
for i in range(1,ny-1):
	for j in range(1,nx-1):
		dist = np.sqrt((x[j]-center[0])**2 + (y[i]-center[1])**2)
		if dist < radius:
			Solid[i,j] = True
			#rho[i,j] = 0

# Iteration
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
			phi_new[i,j] = tools.calc_phi(phi, rho, rho0, i, j, nx, ny, dx, dy, Solid)
	# Last column boundary condition (normal outflow)
	phi_new[:,nx-1] = phi[:,nx-2]
	#print(np.round(phi*100)/100)
	# Meassure error
	error = np.max(np.abs(phi_new - phi))

	for i in range(1, ny-1): # rows
		for j in range (1, nx-1): # columns
			Vx[i,j], Vy[i,j] = tools.calc_vel(phi_new, rho, rho0, i, j, nx, ny, dx, dy, Solid, Vin)
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
	print("Iteration %i: maximum error: %2.4e" %(iter, error))

	# Plots
	#ax2.streamplot(xv, yv, Vx, Vy, density=[0.5, 1])
	#im = ax2.pcolormesh(xv, yv, phi, cmap=cmap)
	#plt.pause(0.05)

print(np.max(np.max(V)))

cmap = plt.get_cmap('PiYG')

fig1 = plt.figure()
ax1 = fig1.gca()
ax1.streamplot(xv, yv, Vx, Vy, density=[0.5, 1])
im = ax1.pcolormesh(xv, yv, V, cmap=cmap)
fig1.colorbar(im, ax=ax1)
circ = Circle(center, radius, fill=False)
ax1.add_patch(circ)
ax1.axis("equal")
plt.title("Velocity (m/s)")

fig2 = plt.figure()
ax2 = fig2.gca()
ax2.streamplot(xv, yv, Vx, Vy, density=[0.5, 1])
im = ax2.pcolormesh(xv, yv, phi, cmap=cmap)
fig2.colorbar(im, ax=ax2)
circ = Circle(center, radius, fill=False)
ax2.add_patch(circ)
ax2.axis("equal")
plt.title("Phi")


fig3 = plt.figure()
ax3 = fig3.gca()
ax3.plot(x,Vx[round(0.1*ny),:],'r-o')
ax3.plot(x,Vx[round(0.25*ny),:],'g-o')
ax3.plot(x,Vx[round(0.5*ny),:],'b-o')

plt.show()
