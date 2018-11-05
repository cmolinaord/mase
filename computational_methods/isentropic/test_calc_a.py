import tools
import tmp


L = 1
H = 1
nx = 5
ny = 5
phi, dx, dy = tools.create_mesh(L, H, nx, ny) # Create a null matrix for phi
print("dx=%1.2f" % dx)
print("dy=%1.2f" % dy)
phi1 = phi + 0
rho0 = 2
rho = phi + rho0
rho[:,0] = 3

Vinput = 5.4

print(phi)
print(phi1)

for i in range(0, ny): # rows
	for j in range (0, nx): # columns
		print("i = %i, j = %i, phi = %1.2f" % (i, j, phi[i,j]))
		phi1[i,j], a,b,c,d = tmp.calc_a(phi, rho, rho0, i, j, nx, ny, dx, dy, Vinput)

print(phi1)
