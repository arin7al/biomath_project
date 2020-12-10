import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math
from mpl_toolkits import mplot3d

def ode3(y, t, k1, k2, k3, c):
    u, v, w = y
    dydt = [u * (1 - k1 * v * u / (1 + u))/c, w, c*w - k2 * v * (1 - k3 * (u ** 2) / (1 + u))]
    return dydt

def calc_ode3(args, u0, v0, w0, ts=10, nt=101):
    y0 = [u0, v0, w0]
    t = np.linspace(0, ts, nt)
    sol = odeint(ode3, y0, t, args)
    return sol

k1 = 1
k2 = -0.02
k3 = 1
c = 1

start_u = 0.0
stop_u = 3.0
start_v = 0.0
stop_v = 3.0
start_w = 0.0
stop_w = 3.0
stop_t = 5
num_u = 5
num_v = 5
num_w = 5
num_t = 10

plt.figure()

u0 = 1
v0 = 1
w0 = 1

args = (k1, k2, k3, c)
sol = calc_ode3(args, u0, v0, w0, stop_t, num_t)

ax = plt.axes(projection='3d')

px = [1.618033988749895]
py = [1.618033988749895]
pz = [0.]
ax.scatter3D(px, py, pz, c=pz, marker= '*')
ox = 0.
oy = 0.
oz = 0.
ax.scatter3D(ox, oy, oz, c=oz, marker= '*')

u0_vec = np.linspace(start_u, stop_u, num_u)
v0_vec = np.linspace(start_v, stop_v, num_v)
w0_vec = np.linspace(start_w, stop_w, num_w)
for u0 in u0_vec:
    for v0 in v0_vec:
        for w0 in w0_vec:
            ax.scatter3D(u0, v0, w0, c=w0, cmap='hsv')
            sol = calc_ode3(args, u0, v0, w0, stop_t, num_t)
            u = sol[:, 0]
            v = sol[:, 1]
            w = sol[:, 2]
            ax.plot3D(u, v, w, 'b', label='system trajectories')
#phase planes

u2 = np.outer(np.linspace(0.1, 3.1, 30), np.ones(30))
u3 = np.outer(np.linspace(0, 3, 30), np.ones(30))
v2 = (1+u2)/u2
v3 = u3.copy().T # transpose
w2 = u2.copy().T
w3 = k2/c * v3*(1 - k3*(u3**2/(1+u3)))

ax.plot_surface(u2, v2, w2,cmap='viridis', edgecolor='none')
ax.plot_surface(u3, v3, w3,cmap='viridis', edgecolor='none')



ax.set_xlim3d(0, stop_u + 1)
ax.set_ylim3d(0,stop_v + 1)
ax.set_zlim3d(0,stop_w + 1)

ax.set_xlabel('u')
ax.set_ylabel('v')
ax.set_zlabel('w')

plt.grid()
plt.show()
