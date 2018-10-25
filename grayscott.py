import matplotlib.pyplot as plt
import numpy as np
import random
from scipy.integrate import odeint


N = 100
Du = 2e-5
Dv = 1e-5
F = 0.025
k = 0.055

EN = 20000

u = np.ones(N)
v = np.zeros(N)
xs = np.arange(N)
xmax = 2.
xmin = 0.
x = np.linspace(xmax, xmin, N)
dx = (xmax - xmin)/N


for i in range(N/4, 3*N/4):
    u[i] = random.random()*0.2 + 0.4
    v[i] = random.random()*0.2 + 0.2


for i in range(EN+1):
    print(i)
    lu = u*0.
    lv = v*0.
    uvv = u*v*v+0.0

    for j in np.arange(1, N-1):
        lu[j] = (u[j+1] - 2.*u[j] + u[j-1])/(dx*dx)
        lv[j] = (v[j+1] - 2.*v[j] + v[j-1])/(dx*dx)
    lu[0] = (u[1] - 2.*u[0] + u[N-1])/(dx*dx)
    lv[0] = (v[1] - 2.*v[0] + v[N-1])/(dx*dx)
    lu[N-1] = (u[0] - 2.*u[N-1] + u[N-2])/(dx*dx)
    lv[N-1] = (v[0] - 2.*v[N-1] + v[N-2])/(dx*dx)

    u += Du*lu - uvv + F*(1. - u)
    v += Dv*lv + uvv - (F + k)*v

    if i % 100 == 0:
        f1 = plt.figure()
        plt.ylim(0, 1)
        plt.plot(xs, u, c='b', label='u')
        plt.plot(xs, v, c='g', label='v')
        plt.title(i)
        plt.legend(["u", "v"])
        nStr = str(i)
        nStr = nStr.rjust(5, '0')
        plt.title(nStr)
        f1.savefig(nStr)
