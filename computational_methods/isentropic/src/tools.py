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

def calc_a(phi, rho, rho0, i, j, nx, ny, dx, dy, Vin, Solid):

	# North border
	if i == 0:
		aN = 0
		phiN = 0
	else:
		aN = ratio_rho(rho0, rho, i, j, "N", dx, dy, Solid) * dx/dy
		phiN = phi[i-1,j]

	# South border
	if i == ny-1:
		aS = 0
		phiS = 0
	else:
		aS = ratio_rho(rho0, rho, i, j, "S", dx, dy, Solid) * dx/dy
		phiS = phi[i+1,j]

	# West border
	if j == 0:
		aW = 0
		phiW = 0
		if i != 0 and i != ny-1:
			aN = Vin
			aS = Vin
	else:
		aW = ratio_rho(rho0, rho, i, j, "W", dx, dy, Solid) * dy/dx
		phiW = phi[i,j-1]

	# East border
	if j == nx-1:
		aE = 0
		phiE = 0
	else:
		aE = ratio_rho(rho0, rho, i, j, "E", dx, dy, Solid) * dy/dx
		phiE = phi[i,j+1]

	aP = aE + aW + aN + aS

	#WARNING
	# Check this
	if aP == 0:
		phiP = 0
	else:
		phiP = (aN*phiN + aS*phiS + aW*phiW + aE*phiE)/aP
	ve = (phiP-phiE)/dy * aE
	vw = (phiW-phiP)/dy * aW
	vn = (phiN-phiP)/dx * aN
	vs = (phiP-phiS)/dx * aS

	return phiP, vn, vs, vw, ve

def ratio_rho(rho0, rho, i, j, dir, dx, dy, Solid):
	if dir == "N":
		rhoX = rho[i-i,j]
		SolidX = Solid[i-1,j]
		d = dy
	elif dir == "S":
		rhoX = rho[i+1,j]
		SolidX = Solid[i+1,j]
		d = dy
	elif dir == "W":
		rhoX = rho[i,j-1]
		SolidX = Solid[i,j-1]
		d = dx
	else:
		rhoX = rho[i,j+1]
		SolidX = Solid[i,j+1]
		d = dx
	SolidP = Solid[i,j]
	rhoP = rho[i,j]
	if SolidP:
		rhoP = 0
	if SolidX:
		rhoX = 0
	if SolidP and SolidX:
		ratio = 0
	else:
		ratio = d / (0.5*d/rho0*rhoP + 0.5*d/rho0*rhoX)
	return ratio
