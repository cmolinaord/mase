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

def calc_phi(phi, rho, rho0, i, j, nx, ny, dx, dy, Vin, Solid):
	print("i=%i, j=%i"%(i,j))
	print("==========")
	aN = ratio_rho(rho0, rho, i, j, "N", dx, dy, Solid) * dx/dy
	aS = ratio_rho(rho0, rho, i, j, "S", dx, dy, Solid) * dx/dy
	aW = ratio_rho(rho0, rho, i, j, "W", dx, dy, Solid) * dy/dx
	aE = ratio_rho(rho0, rho, i, j, "E", dx, dy, Solid) * dy/dx
	phiN = phi[i-1,j]
	phiS = phi[i+1,j]
	phiW = phi[i,j-1]
	phiE = phi[i,j+1]

	aP = aE + aW + aN + aS
	print("a:")
	print("\t%1.1f" % aN)
	print("%1.1f\t%1.1f\t%1.1f" % (aW,aP,aE))
	print("\t%1.1f" % aS)

	phiP = (aN*phiN + aS*phiS + aW*phiW + aE*phiE)/aP

	ve = (phiP-phiE)/dy * aE
	vw = (phiW-phiP)/dy * aW
	vn = (phiN-phiP)/dx * aN
	vs = (phiP-phiS)/dx * aS

	if j == 1:
		phiP = phiW
		vn = Vin
		vs = Vin
	elif j == nx-2:
		phiP = phiE

	print("phiP = %1.1f" % phiP)
	print("v:")
	print("\t%1.1f" % vn)
	print("%1.1f\t\t%1.1f" % (vw,ve))
	print("\t%1.1f" % vs)

	return phiP, vn, vs, vw, ve

def ratio_rho(rho0, rho, i, j, dir, dx, dy, Solid):
	if dir == "N":
		rhoX = rho[i-1,j]
		rhoX0 = rho0[i-1,j]
		SolidX = Solid[i-1,j]
		d = dy
	elif dir == "S":
		rhoX = rho[i+1,j]
		rhoX0 = rho0[i+1,j]
		SolidX = Solid[i+1,j]
		d = dy
	elif dir == "W":
		rhoX = rho[i,j-1]
		rhoX0 = rho0[i,j-1]
		SolidX = Solid[i,j-1]
		d = dx
	else:
		rhoX = rho[i,j+1]
		rhoX0 = rho0[i,j+1]
		SolidX = Solid[i,j+1]
		d = dx
	SolidP = Solid[i,j]
	rhoP = rho[i,j]
	rhoP0 = rho0[i,j]

	print("i=%i,j=%i, dir=%s, rhoX=%1.1f, rhoX0=%1.1f, rhoP=%1.1f, rhoP0=%1.1f, d=%1.1f" % (i,j,dir,rhoX,rhoX0,rhoP,rhoP0,d))

	if SolidX and not SolidP:
		ratio = 2*rhoP0/rhoP
	elif SolidP and not SolidX:
		ratio = 2*rhoX0/rhoX
	elif SolidP and SolidX:
		ratio = 0
	else:
		ratio = d / (0.5*d/rhoP0*rhoP + 0.5*d/rhoX0*rhoX)
	return ratio
