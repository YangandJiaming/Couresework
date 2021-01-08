import numpy as np
np.seterr(invalid='ignore')

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

x1=np.arange(0.43,0.65,0.01)  #range of x1

# relationships of variables
I=((-m*g*np.sin(phi)+k*(x1-d))*((delta-x1)**2)/c)**(1/2)
V=I*R

# get the maximum value of Xe and Ve
print(x1)
print(V)


