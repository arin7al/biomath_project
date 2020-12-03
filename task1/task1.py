import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import math


def ode(y, t, k1, k2, k3):
    u, v = y
    dydt = [u * (1 - k1 * v * u / (1 + u)), k2 * v * (1 - k3 * (u ** 2) / (1 + u))]
    return dydt


def calc_ode(args, u0, v0, ts=10, nt=101):
    y0 = [u0, v0]
    t = np.linspace(0, ts, nt)
    sol = odeint(ode, y0, t, args)
    return sol


def draw_phase_portrait(args, num_u=1, num_v=1, start_u=0., stop_u=5., start_v=0., stop_v=5., stop_t=10, num_t=101):
    u0_vec = np.linspace(start_u, stop_u, num_u)
    v0_vec = np.linspace(start_v, stop_v, num_v)
    for u0 in u0_vec:
        for v0 in v0_vec:
            plt.scatter(u0, v0)
            sol = calc_ode(args, u0, v0, stop_t, num_t)
            u = sol[:, 0]
            v = sol[:, 1]
            trajectories, = plt.plot(u, v, 'b', label='system trajectories')
    # plt.xlabel('u')
    # plt.ylabel('v')
    # plt.grid()
    # plt.show()
    return trajectories


k1 = 1
k2 = -0.02
k3 = 1


def l1(u):
    v = (1 + u) / (u * k1)
    return v


def draw_nullclines(n, start_u=0., stop_u=5., start_v=0., stop_v=5.):
    u = np.linspace(start_u, stop_u, n)
    v = l1(u)
    plt.plot(u, v, 'g', linewidth=3)

    v2 = np.linspace(start_v, stop_v, n)
    u_star = (1 + math.sqrt(1 + 4 * k3)) / (2 * k3)
    u2 = np.repeat(u_star, n)
    isoclines, = plt.plot(u2, v2, 'g', linewidth=3, label='isoclines')
    # plt.plot([1.618], [l1(1.618)], 'r', marker="*", markersize=12)
    return isoclines


start_u = 0.0
stop_u = 4.0
start_v = 0.1
stop_v = 3.0
stop_t = 100
num_u = 10
num_v = 10
num_t = 10000

args = (k1, k2, k3)
start = 0.1
end = 20

isoclines = draw_nullclines(50, start_u, stop_u, start_v, stop_v)
trajectories = draw_phase_portrait(args, num_u, num_v, start_u, stop_u, start_v, stop_v, stop_t, num_t)
plt.plot([0], [0], 'r', marker="*", markersize=12)
fixed_points, = plt.plot([1.618], [l1(1.618)], 'r', marker="*", markersize=12)

plt.xlabel('u')
plt.ylabel('v')
plt.legend([isoclines, trajectories, fixed_points], ['isoclines', 'system trajectories', 'fixed points'])
# plt.xlim([0, stop_u])
# plt.ylim([0, stop_v])
plt.grid()
plt.show()
