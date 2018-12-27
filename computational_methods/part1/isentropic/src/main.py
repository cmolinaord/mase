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
nx = 60
ny = 30

# Load default options
opt = tools.options()

args = sys.argv
if len(args) > 1:
	opt.itermax = int(args[1])
if len(args) > 3:
	nx = int(args[2])
	ny = int(args[3])
if len(args) > 4:
	opt.precission = float(args[4])
elif len(args) > 5:
	print("ERROR: Too many input arguments")
	print("Exiting...")
	exit()

# Domain parameters (m)
L = 20
H = 10
# Solid objects
center = [L/2, H/2]
radius = min(L,H)/6
obs = tools.obstacle(center, radius)

# Printing input data
print("Maximum iterations = %i" % opt.itermax)
print("L = %1.1f" % L)
print("H = %1.1f" % H)
print("[nx,ny] = [%i,%i]" % (nx,ny))
print("Precission = %1.1E" % opt.precission)
print("Compressible? %s" % ('Yes' if opt.compressible else 'No'))
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
w, f = cfd.boundary(w, f, obs)
# Computation of the fluid dynamic
f = cfd.gauss_seidel(w, f, opt)

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
