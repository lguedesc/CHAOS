import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 4000, 10000)
theta = 1
R = 1
L = 1
omega = 1

phi = np.arctan(omega*L/R)

i = (theta/R)*(1/np.sqrt(1 + (omega*t)**2))*np.cos(omega*t - phi)

plt.plot(t,i)
plt.show()