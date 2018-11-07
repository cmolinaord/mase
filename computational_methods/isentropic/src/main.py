import sys
import numpy as np
import tools
import matplotlib.pyplot as plt
from matplotlib.pyplot import Circle

args = sys.argv
ITERMAX = int(args[1])

# Computational parameters
delta = 0.004
nx = 5
ny = 5
iter = 0

# Domain parameters (m)
L = 10
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
Vinput = 5

# Creation of domain
phi, dx, dy = tools.create_mesh(L, H, nx, ny) # Create a null matrix for phi
phi = np.zeros([ny,nx]) + 0.1
phi1 = np.zeros([ny,nx])
p = np.zeros([ny,nx]) + p0
p1 = np.zeros([ny,nx]) + p0
T = np.zeros([ny,nx]) + T0
T1 = np.zeros([ny,nx]) + T0
rho = np.zeros([ny,nx]) + 2.7
rho1 = np.zeros([ny,nx]) + 1.5
Vx = np.zeros([ny,nx])
Vy = np.zeros([ny,nx])
V = np.zeros([ny,nx])
V1 = np.ones([ny,nx])

x = np.linspace(0, L, nx)
y = np.linspace(0, H, ny)
xv, yv = np.meshgrid(x, y)

# Solid objects
center = [L/2, H/2]
radius = min(L,H)/5

Solid, dx, dy = tools.create_mesh(L, H, nx, ny)
Solid[:] = False
#for i in range(0,ny):
#	for j in range(0,nx):
#		dist = np.sqrt((x[j]-center[0])**2 + (y[i]-center[1])**2)
#		if dist < radius:
#			Solid[i,j] = True
#			rho[i,j] = 0

# Boundary conditions
#phi[0, :] = 1

# Iteration
incr = 0.5

while incr > delta and iter < ITERMAX:
	iter += 1
	for i in range(1, ny-1): # rows
		for j in range (0, nx): # columns
			phi1[i,j], vn, vs, vw, ve = tools.calc_a(phi, rho1, rho, i, j, nx, ny, dx, dy, Vinput, Solid)
			Vx[i,j] = 0.5*(vn + vs)
			Vy[i,j] = 0.5*(ve + vw)
			V1[i,j]	= tools.np.sqrt(Vx[i,j]**2 + Vx[i,j]**2)
			# Energy conservation (calculated temperature)
			T1[i,j] = T[i,j] + 0.5*(V[i,j]**2 - V1[i,j]**2)/c_p
			# Isentropic condition (pressure calculated)
			p1[i,j] = p[i,j] * (T1[i,j]/T[i,j])**gamma_exp
	print("rho")
	print(rho)
	rho = rho1
	print("rho = rho1")
	print(rho)
	rho1 = tools.density(p1, T1, R)
	print("rho1 (new)")
	print(rho1)
	incr = np.max(np.abs(phi1 - phi))
	V = V1
	T = T1
	p = p1
	phi = phi1
	print("Iteration %i: maximum difference: %2.4e" %(iter, incr))

cmap = plt.get_cmap('PiYG')

fig1 = plt.figure()
ax1 = fig1.gca()
ax1.streamplot(xv, yv, Vx, Vy, density=[0.5, 1])
im = ax1.pcolormesh(xv, yv, V, cmap=cmap)
fig1.colorbar(im, ax=ax1)
circ = Circle(center, radius)
ax1.add_patch(circ)
ax1.axis("equal")

fig3 = plt.figure()
ax3 = fig3.gca()
ax3.plot(x,Vx[round(0.1*ny),:],'r-o')
ax3.plot(x,Vx[round(0.25*ny),:],'g-o')
ax3.plot(x,Vx[round(0.5*ny),:],'b-o')

plt.show()
