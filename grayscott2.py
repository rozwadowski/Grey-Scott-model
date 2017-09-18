#!/usr/bin/env python
# -*- coding: utf-8 -*-
# P. Rozwadowski
import math
import numpy as np
import random
import matplotlib
import matplotlib.pyplot as plt
from PIL import Image
import scipy as scipy
from scipy import ndimage

#const
N=200
EN = 15000

Du = 2.*1e-5
Dv = 1e-5
F = 0.0220
k= 0.0470

xmax=2.
xmin=0.
x=np.linspace(xmax,xmin,N)
dx=(xmax-xmin)/N

#init
u = np.ones((N,N))
v = np.zeros((N,N))

stencil = np.array([[0, 1, 0],[1, -4, 1], [0, 1, 0]])

for i in range(N/4,3*N/4):
	for j in range(N/4,3*N/4):
		u[i][j] = random.random()*0.2+0.4
		v[i][j] = random.random()*0.2+0.2

for j in range(EN+1):
	print j
	lu=np.zeros((N,N))
	lv=np.zeros((N,N))
	uvv = u*v*v
	# laplasjan
	lu = ndimage.convolve(u, stencil)/(dx**2)
	lv = ndimage.convolve(v, stencil)/(dx**2)
	# laplasjan na brzegu
	for i in np.arange(1,N-1):
		lu[i][0]=(u[i+1][0]+u[i-1][0]-4.*u[i][0]+u[i][N-1]+u[i][1])/(dx*dx)
		lv[i][0]=(v[i+1][0]+v[i-1][0]-4.*v[i][0]+v[i][N-1]+v[i][1])/(dx*dx)
		lu[i][N-1]=(u[i+1][N-1]+u[i-1][N-1]-4.*u[i][N-1]+u[i][N-2]+u[i][0])/(dx*dx)
		lv[i][N-1]=(v[i+1][N-1]+v[i-1][N-1]-4.*v[i][N-1]+v[i][N-2]+v[i][0])/(dx*dx)
		lu[0][i]=(u[0][i+1]+u[0][i-1]-4.*u[0][i]+u[N-1][i]+u[1][i])/(dx*dx)
		lv[0][i]=(v[0][i+1]+v[0][i-1]-4.*v[0][i]+v[N-1][i]+v[1][i])/(dx*dx)
		lu[N-1][i]=(u[N-1][i+1]+u[N-1][i-1]-4.*u[N-1][i]+u[N-2][i]+u[0][i])/(dx*dx)
		lv[N-1][i]=(v[N-1][i+1]+v[N-1][i-1]-4.*v[N-1][i]+v[N-2][i]+v[0][i])/(dx*dx)

	lu[0][0]=(u[1][0]+u[N-1][0]-4.*u[0][0]+u[0][N-1]+u[0][1])/(dx*dx)
	lv[0][0]=(v[1][0]+v[N-1][0]-4.*v[0][0]+v[0][N-1]+v[0][1])/(dx*dx)
	lu[0][N-1]=(u[1][N-1]+u[N-1][N-1]-4.*u[0][N-1]+u[0][N-2]+u[0][0])/(dx*dx)
	lv[0][N-1]=(v[1][N-1]+v[N-1][N-1]-4.*v[0][N-1]+v[0][N-2]+v[0][0])/(dx*dx)
	lu[N-1][0]=(u[N-1][1]+u[N-1][N-1]-4.*u[N-1][0]+u[N-2][0]+u[0][0])/(dx*dx)
	lv[N-1][0]=(v[N-1][1]+v[N-1][N-1]-4.*v[N-1][0]+v[N-2][0]+v[0][0])/(dx*dx)
	lu[N-1][N-1]=(u[0][N-1]+u[N-2][N-1]-4.*u[N-1][N-1]+u[N-1][N-2]+u[N-1][0])/(dx*dx)
	lv[N-1][N-1]=(v[0][N-1]+v[N-2][N-1]-4.*v[N-1][N-1]+v[N-1][N-2]+v[N-1][0])/(dx*dx)

	u+=Du*lu - uvv +  F *(1.-u)
	v+=Dv*lv + uvv - (F+k)*v
	if j%100.==0.:
		fig=plt.figure()
		ax=fig.add_subplot(111)
		cax = ax.imshow(u, interpolation='nearest')
		cax.set_clim(vmin=0, vmax=1)
		cbar = fig.colorbar(cax, ticks=[0,0.3, 0.5,1], orientation='vertical')
		nStr=str(j)
		nStr=nStr.rjust(5,'0')
		plt.title("Grey-Scott model "+nStr+ " (Author: P. Rozwadowski)")
		fig.savefig(nStr)
		plt.clf() 
