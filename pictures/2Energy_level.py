# coding: utf-8

import matplotlib.pyplot as plt
import numpy as np
def chi(x):
  for i in x:
    ans = 1j/(1j*i + 1/2)
    return ans


m = np.arange(-2, 2, 0.01)
y = 1 + chi(m)

plt.plot(m, imag(y), m, real(y))
plt.show()