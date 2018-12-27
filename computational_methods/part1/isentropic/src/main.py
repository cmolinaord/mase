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
import tools
import cfd
import plot
import const as c
import cli

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

tic = time.perf_counter()

# Inlet face
Vin = 3

# Creation of the World
w = tools.world(L, H, nx, ny)
# Print input parameters
cli.input(w, opt)
# Creation of the fluid
f = tools.fluid(w, c.p0, c.T0, Vin)
# Boundary conditions computation
w, f = cfd.boundary(w, f, obs)
# Computation of the fluid dynamic
f = cfd.gauss_seidel(w, f, opt)

toc = time.perf_counter() - tic
print("Elapsed time: %1.2fs" % toc)

plot.velocity(w, f, obs)
plot.density(w, f, obs)
plot.phi(w, f, obs)
plot.show()
