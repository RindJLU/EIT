# coding: utf-8
import numpy as np
import matplotlib.pyplot as plt

# define some constant:
i = 1j
n = pow(10, 19)
e = 1.6*pow(10, -19)
h = 1.05*pow(10, -34)
gama = pow(10, 5)
int_ab = 1.0*pow(10, -30)
eps0 = 8.854*pow(10, -12)
x = np.arange(-2, 2, 0.01)

eps_r =(i*n*pow(int_ab, 2)*gama)/(eps0*h*(1 + i*x))

plt.plot(x, np.imag(eps_r))
plt.plot(x, np.real(eps_r))
plt.show()