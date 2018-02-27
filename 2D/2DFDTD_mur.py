import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import os
import datetime

simulation_time=datetime.datetime.today()
folder_name = "2DFDTD_mur_{0:%Y%m%d-%H%M%S}".format(simulation_time)
os.mkdir(folder_name)

#FDTDの式
def H1_n(i,j):
    if j==step_j-1:
        return H1[i][j]-params[0]*(0-E3[i][j])
    else:
        return H1[i][j]-params[0]*(E3[i][j+1]-E3[i][j])
def H2_n(i,j):
    if i==step_i-1:
        return H2[i][j]+params[0]*(0-E3[i][j])
    else:
        return H2[i][j]+params[0]*(E3[i+1][j]-E3[i][j])
def E3_n(i,j): 
    if j==0:
        return E3[i][j+1]+params[3]*(E3_n(i,j+1)-E3[i][j])
    elif j==step_j-1:
        return past[i][j-1]+params[3]*(E3[i][j-1]-E3[i][j])
    elif i==0:
        return E3[i+1][j]+params[3]*(E3_n(i+1,j)-E3[i][j])
    elif i==step_i-1:
        return past[i-1][j]+params[3]*(E3[i-1][j]-E3[i][j])
    else:
        return params[1]*E3[i][j]+params[2]*(H2[i][j]-H2[i-1][j]-H1[i][j]+H1[i][j-1])
def gausiannpalse(t):
    to=0.5e-9
    a=(4/to)**2
    return math.pow(math.e,-1*a*(t-to)**2)
def E(y,x):
    return E[i][j]
fig = plt.figure()    
#初期化、今回は異物を空気とする
step_i=161
step_j=201
timestep=500
H1=np.zeros((step_i,step_j))
H2=np.zeros((step_i,step_j))
E3=np.zeros((step_i,step_j))
ID=np.zeros((step_i,step_j))
params_=np.zeros((step_i,step_j,3))
params=np.zeros(4)
dx=2.0e-3#f=250,dx=λ*0.25
dt=3.0e-12
ε=66.41*(1e-12)
ε0=8.855*(1e-12)
μ=1.257*(1e-6)
σ=0.000000000001
c=300000000/np.sqrt(ε/ε0)
params[0]=dt/(μ*dx)
params[1]=(1-σ*dt/2/ε)/(1+σ*dt/2/ε)
params[2]=(dt/(dx*ε))/(1+σ*dt/2/ε)
params[3]=(c*dt-dx)/(c*dt+dx)
ims = []

#電磁場の計算
past=[]
for n in range(timestep+1):
    print(n)
    T=n*dt
    E3[80][100]+=gausiannpalse(T)
    past=E3.copy()
    for i in range(step_i):
        for j in range(step_j):
            E3[i][j]=E3_n(i,j)
    
    for i in range(step_i-1):
        for j in range(step_j-1):
            H2[i][j]=H2_n(i,j)
            H1[i][j]=H1_n(i,j)
    if n%20==0:
        
        """
        #静止画表示
        ax = fig.add_subplot(111)
        cax = ax.imshow(E3, interpolation='none',vmin=-0.5,vmax=0.5,animated=True)
        fig.colorbar(cax)
        """
        #動画
        im = plt.imshow(E3,vmin=-0.3,vmax=0.3,animated=True)
        ims.append([im])

ani = animation.ArtistAnimation(fig,ims, interval=50,
                                repeat_delay=1000)

ani.save('{}/2DFDTD_mur.mp4'.format(folder_name), writer="ffmpeg")

#        ims.append([cax])


#        filename = "{}\\output{}ps.png".format(folder_name,n)
#        fig.savefig(filename)
#ani = animation.ArtistAnimation(fig, ims)
#ani.save('anim.mp4', writer="ffmpeg")
 
