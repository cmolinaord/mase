import matplotlib.pyplot as plt
from matplotlib.pyplot import Circle

cmap = plt.get_cmap('jet')

def velocity(w, f, o):
	fig1 = plt.figure()
	ax1 = fig1.gca()
	ax1.streamplot(w.xv, w.yv, f.Vx, f.Vy, density=[0.5, 1])
	im = ax1.pcolormesh(w.xv, w.yv, f.V, cmap=cmap)
	fig1.colorbar(im, ax=ax1)
	circ = Circle(o.c, o.r, fill=False)
	ax1.add_patch(circ)
	ax1.axis("equal")
	plt.title("Velocity (m/s)")

def density(w, f, o):
	fig1 = plt.figure()
	ax1 = fig1.gca()
	ax1.streamplot(w.xv, w.yv, f.Vx, f.Vy, density=[0.2, 0.5])
	im = ax1.pcolormesh(w.xv, w.yv, f.rho, cmap=cmap)
	fig1.colorbar(im, ax=ax1)
	circ = Circle(o.c, o.r, fill=False)
	ax1.add_patch(circ)
	ax1.axis("equal")
	plt.title("Density (Kg m^-3)")

def temperature(w, f, o):
	fig1 = plt.figure()
	ax1 = fig1.gca()
	ax1.streamplot(w.xv, w.yv, f.Vx, f.Vy, density=[0.2, 0.5])
	im = ax1.pcolormesh(w.xv, w.yv, f.T, cmap=cmap)
	fig1.colorbar(im, ax=ax1)
	circ = Circle(o.c, o.r, fill=False)
	ax1.add_patch(circ)
	ax1.axis("equal")
	plt.title("Temperature (K)")

def pressure(w, f, o):
	fig1 = plt.figure()
	ax1 = fig1.gca()
	ax1.streamplot(w.xv, w.yv, f.Vx, f.Vy, density=[0.2, 0.5])
	im = ax1.pcolormesh(w.xv, w.yv, f.p, cmap=cmap)
	fig1.colorbar(im, ax=ax1)
	circ = Circle(o.c, o.r, fill=False)
	ax1.add_patch(circ)
	ax1.axis("equal")
	plt.title("Pressure (Pa)")

def phi(w, f, o):
	fig1 = plt.figure()
	ax1 = fig1.gca()
	im = ax1.contourf(w.xv, w.yv, f.phi, levels=10)
	fig1.colorbar(im, ax=ax1)
	circ = Circle(o.c, o.r, fill=False)
	ax1.add_patch(circ)
	ax1.axis("equal")
	plt.title("Phi")

def show():
	plt.show()
