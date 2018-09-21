
	# - - - - - - - - - - - - - - - - - -
	#                                    |
	#     .····.                         |
	#    ·      ·                        |
	#    ·      ·                        |
	#     '·..·'                         |
	#                                    |
	# # # # # # # # # # # # # # # # # # #

# # = Adiabathic
# | = non-adiabatic

import numpy as np
import matplotlib.pyplot as plt

# Physical parameters
# ******************
# Domain
Lx = 20.0
Ly = 15.0
# Cylinder
Cc = [7,5.6] # Center coordinates (x,y)
Cr = 3.4 # Radius of Cylinder

# Material parameters
# thermal conductivity of materials in W*m^-1*K^-1
kM = 1.0
kC = 6.0

# External thermal conductivity
alpha_ext =1
T_ext = 15.0
# Cylinder heat generation
gV = 2.0

# Mesh parameters
# ******************
Nx = 200
Ny = 150

# Precision of final temeprature (values in K)
delta = 0.1
# Initial temperature for iterations
T0 = 10.0

# Start of mesh calculations
dx = Lx / Nx
dy = Ly / Ny

XX = np.arange(0.5*dx,Lx,dx)
XX = np.append([0],XX)
XX = np.append(XX,[Lx])

YY = np.arange(0.5*dy,Ly,dy)
YY = np.append([0],YY)
YY = np.append(YY,[Ly])

x, y = np.meshgrid(XX, YY)

# Matrices
T = T0 * np.ones([Ny+2, Nx+2])
k = kM * np.ones([Ny+2, Nx+2])
g = np.zeros([Ny+2, Nx+2])

# Inside the Cylinder
for i in range(0, Ny + 1): # row
	for j in range(0, Nx + 1): # column
		d = np.abs(np.sqrt(np.power((x[i,j] - Cc[0]),2) + np.power((y[i,j] - Cc[1]),2)))
		if d < Cr: # If inside the Cylinder
			k[i,j] = kC
			g[i,j] = gV

			



















#plt.pcolor(x,y,g)
#plt.colorbar()
#plt.show()
