import sys
import numpy as np
import tools
import matplotlib.pyplot as plt
from matplotlib.pyplot import Circle

args = sys.argv
ITERMAX = int(args[1])

# Computational parameters
delta = 0.005
nx = 200
ny = 200
iter = 0

# Domain parameters (m)
L = 10
H = 10

# Physical constants
R = 287.058
c_p = 0.0010035 # Isobaric mass heat capacity (c_p) (J*kg^-1*K^-1)
c_p = 1.0035
gamma = 1.4 # Heat capacity ratio for dry air at room temperature

# Physical input values
T0 = 273 # Initial temperature (Kelvin)
p0 = 101325 # Initial pressure (Pa)
rho0 = tools.density(p0, T0, R) # Intial density (kg/m^3)

# Inlet face
Vinput = 10

# Creation of domain
phi, dx, dy = tools.create_mesh(L, H, nx, ny) # Create a null matrix for phi
phi = np.zeros([nx,ny])
phi1 = np.ones([nx,ny])
p = np.ones([nx,ny]) + p0
T = np.ones([nx,ny]) + T0
rho = np.ones([nx,ny]) * p0
Vx = np.ones([nx,ny])
Vy = np.zeros([nx,ny])
V = np.ones([nx,ny])

x = np.linspace(0, L, nx)
y = np.linspace(0, H, ny)
xv, yv = np.meshgrid(x, y)

# Solid objects
center = [L/2, H/2]
radius = min(L,H)/5

Solid, dx, dy = tools.create_mesh(L, H, nx, ny)
Solid[:] = False
for i in range(0,ny):
	for j in range(0,nx):
		dist = np.sqrt((x[j]-center[0])**2 + (y[i]-center[1])**2)
		if dist < radius:
			Solid[i,j] = True
			rho[i,j] = 0

# Boundary conditions
phi[0, :] = 1

# Iteration
incr = 0.5

while incr > delta and iter < ITERMAX:
	iter += 1
	for i in range(0, ny): # rows
		for j in range (0, nx): # columns
			phi1[i,j], vn, vs, vw, ve = tools.calc_a(phi, rho, rho0, i, j, nx, ny, dx, dy, Vinput, Solid)
			Vx[i,j] = 0.5*(vn + vs)
			Vy[i,j] = 0.5*(ve + vw)
			V[i,j]	= tools.np.sqrt(Vx[i,j]**2 + Vx[i,j]**2)
			# Energy conservation (calculated temperature)
			T[i,j] = T0 + 0.5*(Vinput**2 - V[i,j]**2)/c_p
			# Isentropic condition (pressure calculated)
			p[i,j] = p0 * (T[i,j]/T0)**(gamma/(gamma-1))
			rho[i,j] = tools.density(p[i,j], T[i,j], R)
	#phi1[:,0] = 1
	incr = np.max(np.abs(phi1 - phi))
	phi[:,:] = phi1[:,:]
	print("Iteration %i: maximum difference: %2.4e" %(iter, incr))


fig1 = plt.figure()
ax1 = fig1.gca()
ax1.streamplot(xv, yv, Vx, Vy, density=[0.5, 1])
circ = Circle(center, radius)
ax1.add_patch(circ)
ax1.axis("equal")

cmap = plt.get_cmap('PiYG')

fig2 = plt.figure()
ax2 = fig2.gca()
im = ax2.pcolormesh(xv, yv, V, cmap=cmap)
fig2.colorbar(im, ax=ax2)

fig3 = plt.figure()
ax3 = fig3.gca()
ax3.plot(x,Vx[round(0.1*ny),:],'r-o')
ax3.plot(x,Vx[round(0.25*ny),:],'g-o')
ax3.plot(x,Vx[round(0.5*ny),:],'b-o')

plt.show()
