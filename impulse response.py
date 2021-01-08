import numpy as np
import sympy as sym


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
s, t = sym.symbols('s, t')
transfer_function = 1 /(O*(s**3) + P*(s**2) + J*s + K)

# Impulse response and Step response of the system
impulse_response=sym.inverse_laplace_transform(transfer_function * 1 , s, t)
step_response=sym.inverse_laplace_transform(transfer_function * (1/s) , s, t)

sym.pprint(impulse_response.simplify())
sym.pprint(step_response.simplify())

