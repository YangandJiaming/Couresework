import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import sympy as sym

m = 0.425
g = 9.81
d = 0.42
delta = 0.65
r = 0.125
R = 53
L0 = 0.12
L1 = 0.25
alpha = 1.2
c = 6.815
k = 1880
b = 10.4
phi = 5.25
V = 39.7


def ballmotion(_t, z):
    x1 = z[0]
    x2 = z[1]
    I = z[2]
    return [x2,
            (5*m/7)*(m*g*np.sin(5.25)-k*(x1-d)-b*x2+c*((I)**2)/((delta-x1)**2)),
            (V-I*R)/(L0+L1*np.exp(-alpha*(delta-x1)))]


z_initial = [0.4876, 0, 0.75]
t_final = 2
num_points = 1000
sol = solve_ivp(ballmotion,
                [0, t_final],
                z_initial,
                t_eval=np.linspace(0, t_final, num_points))

x1_simulation = sol.y[0]

plt.plot(sol.t, sol.y[0].T)
plt.show()