import numpy as np

def create_mesh(L, H, nx, ny):
	A = np.zeros([ny, nx])
	dx = L/nx
	dy = H/ny
	return A, dx, dy

def density(p, T, R):
	# Calculate density of the gas given p (pressure), T (temperature) and
	# the specific gas constant R
	rho = p / T / R
	return rho

def calc_phi(phi, rho, rho0, i, j, nx, ny, dx, dy, Solid):

	aN = ratio_rho(rho, rho0, i, j, "N", dx, dy, Solid) * dx/dy
	aS = ratio_rho(rho, rho0, i, j, "S", dx, dy, Solid) * dx/dy
	aW = ratio_rho(rho, rho0, i, j, "W", dx, dy, Solid) * dy/dx
	aE = ratio_rho(rho, rho0, i, j, "E", dx, dy, Solid) * dy/dx
	phiN = phi[i+1,j]
	phiS = phi[i-1,j]
	phiW = phi[i,j-1]
	phiE = phi[i,j+1]

	aP = aE + aW + aN + aS
	# If fully surrounded by Solid
	phi_solid = 0.5
	if Solid[i,j] == True:
		phiP = phi[int(phi_solid*ny),0]
	else:
		phiP = (aN*phiN + aS*phiS + aW*phiW + aE*phiE)/aP

	return phiP

def calc_vel(phi, rho, rho0, i, j, nx, ny, dx, dy, Solid, Vin):
	phiP = phi[i,j]
	phiN = phi[i+1,j]
	phiS = phi[i-1,j]
	phiW = phi[i,j-1]
	phiE = phi[i,j+1]

	vn = (phiN-phiP)/dy * ratio_rho(rho, rho0, i, j, "N", dx, dy, Solid)
	vs = (phiP-phiS)/dy * ratio_rho(rho, rho0, i, j, "S", dx, dy, Solid)
	vw = (phiW-phiP)/dx * ratio_rho(rho, rho0, i, j, "W", dx, dy, Solid)
	ve = (phiP-phiE)/dx * ratio_rho(rho, rho0, i, j, "E", dx, dy, Solid)

	vx = 0.5*(vn + vs)
	vy = 0.5*(ve + vw)
	return vx, vy

def ratio_rho(rho, rho0, i, j, dir, dx, dy, Solid):
	# Central node
	SolidP = Solid[i,j]
	rhoP = rho[i,j]
	# Neighbor node
	if dir == "N":
		rhoX = rho[i+1,j]
		SolidX = Solid[i+1,j]
		d = dy
	elif dir == "S":
		rhoX = rho[i-1,j]
		SolidX = Solid[i-1,j]
		d = dy
	elif dir == "W":
		rhoX = rho[i,j-1]
		SolidX = Solid[i,j-1]
		d = dx
	elif dir == "E":
		rhoX = rho[i,j+1]
		SolidX = Solid[i,j+1]
		d = dx

	if SolidX and not SolidP:
		ratio = 2*rho0/rhoP
	elif SolidP and not SolidX:
		ratio = 2*rho0/rhoX
	elif SolidP and SolidX:
		ratio = 0
	else:
		ratio = d / (0.5*d/rho0*rhoP + 0.5*d/rho0*rhoX)
	return ratio
