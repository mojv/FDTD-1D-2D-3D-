# -*- coding: utf-8 -*-
"""
FDTD1D　
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import os
import datetime
import math

simulation_time=datetime.datetime.today()
folder_name = "1DFDTD_noABC_{0:%Y%m%d-%H%M%S}".format(simulation_time)
import matplotlib.animation as animation
os.mkdir(folder_name)

def E(i):
    if i==0:
        return params[0]*Ex[i]-params[1]*(Hy[i]-0)
    elif i==6000:
        return params[0]*Ex[i]-params[1]*(0-Hy[i-1])
    else:
        return params[0]*Ex[i]-params[1]*(Hy[i]-Hy[i-1])
    
def H(i):
    return Hy[i]-params[2]*(Ex[i+1]-Ex[i])
    
def gausiannpalse(t):
    to=2e-9
    a=(4/to)**2
    return math.pow(math.e,-1*a*(t-to)**2)

start=time.clock()
Ex=np.zeros(6001)
Hy=np.zeros(6000)
dt=3*(1e-12)
dz=0.001
ε=8.85*(1e-12)
μ=1.257*(1e-6)
σ=0.001
X=np.arange(6001)
params=np.zeros(3)
params[0]=(2*ε-dt*σ)/(2*ε+dt*σ)
params[1]=(2*dt)/(dz*(2*ε+dt*σ))
params[2]=dt/(dz*μ)

sousingen=3000

fig=plt.figure()
ims=[]
for n in range(1,12000):
    t=n*3*(1e-12)
    ti=(n-1)*dt
    Ex[sousingen]+=(gausiannpalse((t+ti)/2))
#    Ex[sousingen]+=gausiannpalse(t)
    for i in range(6001):
        Ex[i]=E(i)
    for i in range(6000):
        Hy[i]=H(i)
    
    
    if n%100==0:
        print(n)
        xmin,xmax=0,6002
        ymin = -1
        ymax = 1
        plt.axis([xmin, xmax, ymin, ymax])
#        plt.plot(X,Ex)
#        plt.show()
        line,=plt.plot(X,Ex,"r")
        ims.append([line])

ani = animation.ArtistAnimation(fig, ims)
    
ani.save('{}/1D_noABC.mp4'.format(folder_name), writer="ffmpeg")

#        filename = "{}\\output{}ps.png".format(folder_name,2*n)
#        plt.savefig(filename)


end=time.clock()
print(end-start)