import numpy as np
import matplotlib.pyplot as plt

# Reference
V = 1

# Parameters 
Lp = 0.3
g = 9.81
mp = 0.3
ms = 0.6
kz = 300
kx = 60
Lem = 1e-2
Rem = 10
theta_em = 9e-4

omega_z = np.sqrt(kz/ms)

# Normalized
chi_em = (theta_em*V)/(mp*Lp*Lp*omega_z*omega_z)
varphi_em = Rem/(Lem*omega_z)
kappa_em = theta_em*Rem/(Lem*V)

print(f"chi_em = {chi_em}")
print(f"varphi_em = {varphi_em}")
print(f"kappa_em = {kappa_em}")


plt.show()