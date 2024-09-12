#
# Sensitivity of the Lorenz system to the coefficients
#
# Author: Dr Chennakesava Kadapa
#
# Date: 12-Sep-2024
#
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.cm import rainbow
import numpy as np
from scipy.integrate import odeint, solve_ivp
from scipy.io import loadmat


def Lorenz1(x, t):
    return [
        -10.0*x[0] + 10.0*x[1],
         28.0*x[0] - x[1] - x[0]*x[2],
         x[0]*x[1] - 8.0*x[2]/3.0
    ]

def Lorenz2(x, t):
    return [
        -10.0*x[0] + 10.0*x[1],
         28.0*x[0] - x[1] - x[0]*x[2],
         x[0]*x[1] - 2.667*x[2]
    ]

dt = 0.01
t = np.arange(0, 20, dt)

x0 = [12.345678910111213,
      13.121110987654321,
      22.345678910203040]

sol_L1 = odeint(Lorenz1, x0, t)
sol_L2 = odeint(Lorenz2, x0, t)



fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
axs[0].plot(t, sol_L1[:, 0], "b",   label=r"$\beta=8.0/3$", lw=2, alpha=0.5)
axs[0].plot(t, sol_L2[:, 0], "k--", label=r"$\beta=2.667$", lw=2)
axs[0].legend()
axs[0].set_ylabel("$x$", fontsize=14)
axs[0].set_xlim([0.0,20.0])

axs[1].plot(t, sol_L1[:, 1], "b",   label=r"$\beta=28/3$", lw=2, alpha=0.5)
axs[1].plot(t, sol_L2[:, 1], "k--", label=r"$\beta=2.667$", lw=2)
axs[1].set_ylabel("$y$", fontsize=14)

axs[2].plot(t, sol_L1[:, 2], "b",   label=r"$\beta=28/3$", lw=2, alpha=0.5)
axs[2].plot(t, sol_L2[:, 2], "k--", label=r"$\beta=2.667$", lw=2)
#axs[2].legend()
axs[2].set_xlabel("t", fontsize=14)
axs[2].set_ylabel("$z$", fontsize=14)

plt.savefig("plot.png")















