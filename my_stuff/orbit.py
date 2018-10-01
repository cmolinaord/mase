import numpy as np
from astropy import units as u
from astropy import constants as const
from matplotlib import pyplot as plt

def norm(vec):
	n = np.sqrt(vec[0]**2 + vec[1]**2)
	return n

def accel(X):
	r3 = norm(X.to(u.m))**3
	a = - X / r3 * const.GM_earth
	return a.to(u.m/u.s**2)

R_e = const.R_earth.to(u.km).value

# Earth fixed position
Earth = np.array([0, 0])
circle = plt.Circle((0,0), R_e, color='blue')

# Time to integrate in seconds
t_end = 2 * u.hour
dt = 1 * u.s
N = int(t_end.to(dt.unit) / dt)
t = np.linspace(0, t_end.to(u.second), N)
print("Total integration time = %i%s" % (t_end.value,t_end.unit))
print("dt = %i%s" % (dt.value,dt.unit))
print("Number of itereations = %i" % N)

# Initial height over Earth surface
h0 = 400 * u.km
# Initial position
x0 = h0 + const.R_earth
y0 = 0
# Initial velocity module
v0 = 8000 * u.m/u.s

X = np.zeros([N, 2]) * u.km
X[0][0] = x0
X[0][1] = 0

V = np.zeros([N, 2]) * u.m / u.s
V[0][0] = 0
V[0][1] = v0

a = np.zeros([N, 2]) * u.m / u.s**2
a[0] = accel(X[0])

for i in range(1, N):
	dt = t[i] - t[i-1]
	X[i] = X[i-1] + V[i-1] * dt + 0.5 * a[i-1] * dt**2
	V[i] = V[i-1] + a[i-1] * dt
	a[i] = accel(X[i-1])

plt.plot(X[:,0],X[:,1])
fig = plt.gcf()
ax = fig.gca()
ax.add_artist(circle)
ax.set_aspect('equal', 'box')
plt.show()
