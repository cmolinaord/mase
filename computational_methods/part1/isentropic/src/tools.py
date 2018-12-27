import numpy as np
import const as c

class world(object):
	def __init__(self, L, H, nx, ny):
		self.L = L
		self.H = H
		self.nx = nx
		self.ny = ny
		self.dx = L/nx
		self.dy = H/ny
		self.x = np.linspace(0, L, nx)
		self.y = np.linspace(0, H, ny)
		self.xv, self.yv = np.meshgrid(self.x, self.y)
		# Creation of solid
		self.solid = self.create_matrix(0)
	def create_matrix(self, init):
		A = np.zeros([self.ny, self.nx]) + init
		return A

class fluid(object):
	def __init__(self, w, p0, T0, Vin):
		self.phi 	= w.create_matrix(0)
		self.phi_1 	= w.create_matrix(0)
		self.p 	= w.create_matrix(p0)
		self.T	= w.create_matrix(T0)
		self.rho	= w.create_matrix(density(p0, T0, c.R))
		self.Vin	= Vin
		self.Vx	= w.create_matrix(Vin)
		self.Vy	= w.create_matrix(0)
		self.V	= w.create_matrix(0)

class obstacle(object):
	def __init__(self, c, r):
		self.c = c
		self.r = r

class options(object):
	precission = 1e-3
	itermax = 200
	compressible = True
	verbose = True

def density(p, T, R):
	# Calculate density of the gas given p (pressure), T (temperature) and
	# the specific gas constant R
	rho = p / T / R
	return rho

def calc_phi(phi, rho, rho0, i, j, world):

	aN = ratio_rho(rho, rho0, i, j, "N", world) * world.dx/world.dy
	aS = ratio_rho(rho, rho0, i, j, "S", world) * world.dx/world.dy
	aW = ratio_rho(rho, rho0, i, j, "W", world) * world.dy/world.dx
	aE = ratio_rho(rho, rho0, i, j, "E", world) * world.dy/world.dx
	phiN = phi[i+1,j]
	phiS = phi[i-1,j]
	phiW = phi[i,j-1]
	phiE = phi[i,j+1]

	aP = aE + aW + aN + aS
	# If fully surrounded by solid
	phi_solid = 0.5
	if world.solid[i,j] == True:
		phiP = phi[int(phi_solid*world.ny),0]
	else:
		phiP = (aN*phiN + aS*phiS + aW*phiW + aE*phiE)/aP

	return phiP

def calc_vel(phi, rho, rho0, i, j, world, Vin):
	phiP = phi[i,j]
	phiN = phi[i+1,j]
	phiS = phi[i-1,j]
	phiW = phi[i,j-1]
	phiE = phi[i,j+1]

	vn = (phiN-phiP)/world.dy * ratio_rho(rho, rho0, i, j, "N", world)
	vs = (phiP-phiS)/world.dy * ratio_rho(rho, rho0, i, j, "S", world)
	vw = (phiW-phiP)/world.dx * ratio_rho(rho, rho0, i, j, "W", world)
	ve = (phiP-phiE)/world.dx * ratio_rho(rho, rho0, i, j, "E", world)

	vx = 0.5*(vn + vs)
	vy = 0.5*(ve + vw)
	return vx, vy

def ratio_rho(rho, rho0, i, j, dir, world):
	# Central node
	solidP = world.solid[i,j]
	rhoP = rho[i,j]
	# Neighbor node
	if dir == "N":
		rhoX = rho[i+1,j]
		solidX = world.solid[i+1,j]
		d = world.dy
	elif dir == "S":
		rhoX = rho[i-1,j]
		solidX = world.solid[i-1,j]
		d = world.dy
	elif dir == "W":
		rhoX = rho[i,j-1]
		solidX = world.solid[i,j-1]
		d = world.dx
	elif dir == "E":
		rhoX = rho[i,j+1]
		solidX = world.solid[i,j+1]
		d = world.dx

	if solidX and not solidP:
		ratio = 2*rho0/rhoP
	elif solidP and not solidX:
		ratio = 2*rho0/rhoX
	elif solidP and solidX:
		ratio = 0
	else:
		ratio = d / (0.5*d/rho0*rhoP + 0.5*d/rho0*rhoX)
	return ratio
