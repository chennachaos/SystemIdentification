#
# Sensitivity of the Lorenz system to the initial conditions
#
# Author: Dr Chennakesava Kadapa
#
# Date: 12-Sep-2024
#
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.integrate import odeint

def Lorenz(x, t):
    return [
        -10.0*x[0] + 10.0*x[1],
         28.0*x[0] - x[1] - x[0]*x[2],
         x[0]*x[1] - 8.0*x[2]/3.0
    ]

dt = 0.01
t = np.arange(0, 20, dt)

x0_IC1 = [12.345678910111213,
          13.121110987654321,
          22.345678910203040]

sol_IC1 = odeint(Lorenz, x0_IC1, t)

x0_IC2 = [12.346678910111213,
          13.121110987654321,
          22.345678910203040]
sol_IC2 = odeint(Lorenz, x0_IC2, t)

fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True)
axs[0].plot(t, sol_IC1[:, 0], "b",   label="IC1", lw=2, alpha=0.5)
axs[0].plot(t, sol_IC2[:, 0], "k--", label="IC2", lw=2)
axs[0].legend()
axs[0].set_ylabel("$x$", fontsize=14)
axs[0].set_xlim([0.0,20.0])

axs[1].plot(t, sol_IC1[:, 1], "b",   label="IC1", lw=2, alpha=0.5)
axs[1].plot(t, sol_IC2[:, 1], "k--", label="IC2", lw=2)
axs[1].set_ylabel("$y$", fontsize=14)

axs[2].plot(t, sol_IC1[:, 2], "b",   label="IC1", lw=2, alpha=0.5)
axs[2].plot(t, sol_IC2[:, 2], "k--", label="IC2", lw=2)
#axs[2].legend()
axs[2].set_xlabel("t", fontsize=14)
axs[2].set_ylabel("$z$", fontsize=14)

plt.savefig("plot.png")

