#
# Standard regression technique for identifying
# Lorenz system from time series data
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
from scipy import linalg
from copy import deepcopy

# Generate training data

def Lorenz(x, t):
    return [
        -10.0*x[0] + 10.0*x[1],
        28.0*x[0] - x[1] - x[0]*x[2],
        x[0]*x[1] - 8.0*x[2]/3.0
    ]

dt = 0.1
t_train = np.arange(0, 40, dt)
x0_train = [ -8, 8, 27]

sol_train = odeint(Lorenz, x0_train, t_train)
# get data for x, y, z variables
x = sol_train[:,0]
y = sol_train[:,1]
z = sol_train[:,2]

# find the derivative terms
xd = -10.0*x + 10.0*y
yd = 28.0*x - y - x*z
zd = x*y - 8.0*z/3.0

# Fit the model
# Construct matrix A
N = np.shape(t_train)[0]
A = np.zeros( (3*N,30), dtype="float64")
vvec = np.zeros( (3*N,1), dtype="float64")

#for i in range(2,N-2):
for i in range(N):
    ind1 = 3*i
    ind2 = ind1+1
    ind3 = ind1+2

    A[ind1, 0]  = 1.0
    A[ind1, 1]  = x[i]
    A[ind1, 2]  = y[i]
    A[ind1, 3]  = z[i]
    A[ind1, 4]  = x[i]*x[i]
    A[ind1, 5]  = y[i]*y[i]
    A[ind1, 6]  = z[i]*z[i]
    A[ind1, 7]  = x[i]*y[i]
    A[ind1, 8]  = y[i]*z[i]
    A[ind1, 9]  = z[i]*x[i]

    A[ind2, 10:20] = A[ind1, 0:10]
    A[ind3, 20:30] = A[ind1, 0:10]

    vvec[ind1] = xd[i]
    vvec[ind2] = yd[i]
    vvec[ind3] = zd[i]

coeffs = linalg.inv(A.transpose().dot(A)).dot(A.transpose().dot(vvec))

print("Coefficients of [1 x y z x^2 y^2 z^2 xy yz zx] terms")
coeffs = np.reshape(coeffs, (3,10) )
print(coeffs)


def Lorenz_sim(x, t):
    return [
        coeffs[0,0] + coeffs[0,1]*x[0] + coeffs[0,2]*x[1] + coeffs[0,3]*x[2] + coeffs[0,4]*x[0]*x[0]+ coeffs[0,5]*x[1]*x[1]+ coeffs[0,6]*x[2]*x[2] + coeffs[0,7]*x[0]*x[1]+ coeffs[0,8]*x[1]*x[2]+ coeffs[0,9]*x[0]*x[2],
        coeffs[1,0] + coeffs[1,1]*x[0] + coeffs[1,2]*x[1] + coeffs[1,3]*x[2] + coeffs[1,4]*x[0]*x[0]+ coeffs[1,5]*x[1]*x[1]+ coeffs[1,6]*x[2]*x[2] + coeffs[1,7]*x[0]*x[1]+ coeffs[1,8]*x[1]*x[2]+ coeffs[1,9]*x[0]*x[2],
        coeffs[2,0] + coeffs[2,1]*x[0] + coeffs[2,2]*x[1] + coeffs[2,3]*x[2] + coeffs[2,4]*x[0]*x[0]+ coeffs[2,5]*x[1]*x[1]+ coeffs[2,6]*x[2]*x[2] + coeffs[2,7]*x[0]*x[1]+ coeffs[2,8]*x[1]*x[2]+ coeffs[2,9]*x[0]*x[2]
    ]


# Simulate and plot the results

t_sim = np.arange(0, 40, dt)
x_sim = odeint(Lorenz_sim, x0_train, t_sim)

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
axs[0].set_title("With Standard Regression")

axs[1].plot(t_exact, x_exact[:, 1], "b", label="Truth", lw=2, alpha=0.5)
axs[1].plot(t_sim,   x_sim[:, 1], "k--", lw=2)
axs[1].set_ylabel("$y$", fontsize=14)

axs[2].plot(t_exact, x_exact[:, 2], "b", label="Truth", lw=2, alpha=0.5)
axs[2].plot(t_sim,   x_sim[:, 2], "k--", lw=2)
axs[2].set_ylabel("$z$", fontsize=14)

plt.savefig("plot_StdRegr.png")

