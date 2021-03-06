import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


# Given values
m = 0.425
g = 9.81
d = 0.42
delta = 0.65
r = 0.125
R = 53
L0 = 0.12
L1 = 0.025
alpha = 1.2
c = 6.815
k = 1880
b = 10.4
phi = 42.*np.pi/180
x1e=0.478

# Equilibrium points
Ie=(((-m*g*np.sin(phi)+k*(x1e-d))*((delta-x1e)**2)/c))**(1/2)
Ve=Ie*R


# linearised dynamical equations of the system
def ballmotion(_t, z):
    x1 = z[0]
    x2 = z[1]
    I = z[2]
    return [x2,
            -(5*k/7*m)*x1-(5*b/7*m)*x2+(5*c/7*m)*(2*I*Ie/((delta-x1)**2) +2*x1*(Ie ** 2)/((delta-x1) ** 3)),
            -R*I*(L0+L1 * np.exp(-alpha*(delta-x1)))-x1*((Ve-Ie*R)*L1*np.exp(-alpha * (delta - x1)))/((L0+L1*np.exp(-alpha*(delta-x1))))**2]

# Initial value theorem
z_initial = [x1e, 0, Ie]
t_final = 5
num_points = 1000
sol = solve_ivp(ballmotion,
                [0, t_final],
                z_initial,
                t_eval=np.linspace(0, t_final, num_points))
time = sol.t
x1_trajectory = sol.y[0]

# Plot the Trajectory of ball
plt . plot(time, x1_trajectory.T)
plt.xlabel('Time (s)')
plt.ylabel('Ball position (m)')
plt.title('Ball trajectory in linearised system')
plt.grid()
plt . show()
