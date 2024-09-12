#
# SINDy for identifying Lorenz system from time series data
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
from copy import deepcopy

import pysindy as ps

# Generate training data

def lorenz(x, t):
    return [
        -10.0*x[0] + 10.0*x[1],
        28.0*x[0] - x[1] - x[0]*x[2],
        x[0]*x[1] - 8.0*x[2]/3.0
    ]

dt = 0.1
t_train = np.arange(0, 40, dt)
x0_train = [ -8, 8, 27]

sol_train = odeint(lorenz, x0_train, t_train)

# get data for x, y, z variables
x = sol_train[:,0]
y = sol_train[:,1]
z = sol_train[:,2]

# find the derivative terms
xd = -10.0*x + 10.0*y
yd = 28.0*x - y - x*z
zd = x*y - 8.0*z/3.0

x_dot_measured = np.array(
    [lorenz(sol_train[i],0) for i in range(t_train.size)]
)

# Fit the model

poly_order = 5
threshold = 0.0001


# with default options
model = ps.SINDy()

#model = ps.SINDy(
#    optimizer=ps.STLSQ(threshold=threshold),
#    feature_library=ps.PolynomialLibrary(degree=poly_order),
#)

#model = ps.SINDy(
#    optimizer=ps.SR3(threshold=threshold),
#    feature_library=ps.PolynomialLibrary(degree=poly_order),
#)

print("With measured derivative")
model.fit(sol_train, t=dt, x_dot=x_dot_measured)
model.print()

# Simulate and plot the results

t_sim = np.arange(0, 40, dt)
x_sim = model.simulate(x0_train, t_sim)

plot_kws = dict(linewidth=2)

t_exact = deepcopy(t_train)
x_exact = deepcopy(sol_train)
# Code in case we want to simulate for other time instants
#t_exact = np.arange(0, 40, dt)
#x_exact = odeint(Lorenz, x0_train, t_exact)

fig, axs = plt.subplots(3, 1, figsize=(10, 6), sharex=True)

axs[0].plot(t_exact, x_exact[:, 0], "b",   label="Truth", lw=2, alpha=0.5)
axs[0].plot(t_sim,   x_sim[:, 0],   "k--", label="model", lw=2)
axs[0].legend(loc='upper left')
axs[0].set_ylabel("$x$", fontsize=14)
axs[0].set_xlim([0.0,40.0])
axs[0].set_title("With SINDy")

axs[1].plot(t_exact, x_exact[:, 1], "b", label="Truth", lw=2, alpha=0.5)
axs[1].plot(t_sim,   x_sim[:, 1], "k--", lw=2)
axs[1].set_ylabel("$y$", fontsize=14)

axs[2].plot(t_exact, x_exact[:, 2], "b", label="Truth", lw=2, alpha=0.5)
axs[2].plot(t_sim,   x_sim[:, 2], "k--", lw=2)
axs[2].set_ylabel("$z$", fontsize=14)

plt.savefig("plot_SINDy.png")

