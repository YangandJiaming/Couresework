import numpy as np
import sympy as sym
import control as ctrl
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


# Values of Parameters
A=(k-(2*c*(Ie)**2)/((delta-x1e)**3))*(((delta-x1e)**2)/(2*c*Ie))
B=b*((delta-x1e)**2)/(2*c*Ie)
C=(7*m*((delta-x1e)**2))/(10*c*Ie)
D=L0+L1*np.exp(-alpha*(delta-x1e))
E=((Ve-Ie*R)*L1*alpha*np.exp(-alpha*(delta-x1e)))/(L0+L1*np.exp(-alpha*(delta-x1e)))
O=C*D
P=B*D+C*R
J=A*D+B*R
K=A*R+E


# Transfer function of the system
transfer_function = ctrl.TransferFunction([1], [O,P,J,K])

# determine impulse response and step response of system
s, t = sym.symbols('s, t')
transfer_function_1 = 1/(O*(s)**3 + P*(s)**2 + J*s + K)
F_input_impulse = 1
F_input_step = 1/s
impulse_response=sym.inverse_laplace_transform(transfer_function_1 * F_input_impulse, s, t)
step_response=sym.inverse_laplace_transform(transfer_function_1 * F_input_step, s, t)

# Print the response 
sym.pprint(impulse_response.simplify())
sym.pprint(step_response.simplify())

# Impulse response and step response of this system
t,y = ctrl. impulse_response(transfer_function,
                             np.linspace(0, 5, 5000))

t1,y1 = ctrl. step_response(transfer_function,
                             np.linspace(0, 5, 5000))


# Plot the response of system
plt.xlabel('Time (s)')
plt.ylabel('Ball position (m)')
plt.title('Impulse Response and Step Response of System')
plt.plot(t, y,label='impulse response')
plt.plot(t1, y1,label='step response')
plt.legend(loc='upper right')
plt.grid()
plt.show()

