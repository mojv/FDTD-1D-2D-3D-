
"""
FDTD 3d
"""
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import math
from matplotlib.colors import ListedColormap
import os


"""
値
"""
nx=30
ny=30
nz=30
nt=1000
frequency=250.0e-6
"""
ε=8.855*(1e-12)
σ=0.000000000001
μ=1.25663706*1e-6
"""
ε=66.41*(1e-12)
ε0=8.855*(1e-12)
μ=1.257*(1e-6)
σ=0.000000000001

ds=2.0e-3
dt=3.0e-12
c=300000000/np.sqrt(ε/ε0)

domain=np.zeros((nx,ny,nz))
source=domain[int(nx/2)][int(ny/2)][int(nz/2)]
Ex=np.zeros((nx,ny,nz))
Ey=np.zeros((nx,ny,nz))
Ez=np.zeros((nx,ny,nz))
Hx=np.zeros((nx,ny,nz))
Hy=np.zeros((nx,ny,nz))
Hz=np.zeros((nx,ny,nz))

CE=(1-σ*dt/(2*ε))/(1+(σ*dt)/(2*ε))
CEL=(dt/(ds*ε))/(1+σ*dt/2/ε)
CHL=dt/μ/ds
mur=(c*dt-ds)/(c*dt+ds)
def Ex_n(i,j,k):
    if i==0:
        return Ex[i+1][j][k]+mur*(Ex_n(i+1,j,k)-Ex[i][j][k])
    elif i==ny-1:
        return cpEx[i-1][j][k]+mur*(Ex[i-1][j][k]-Ex[i][j][k])
    
    elif j==0:
        return Ex[i][j+1][k]+mur*(Ex_n(i,j+1,k)-Ex[i][j][k])
    elif j==ny-1:
        return cpEx[i][j-1][k]+mur*(Ex[i][j-1][k]-Ex[i][j][k])
    elif k==0:
        return Ex[i][j][k+1]+mur*(Ex_n(i,j,k+1)-Ex[i][j][k])
    elif k==nz-1:
        return cpEx[i][j][k-1]+mur*(Ex[i][j][k-1]-Ex[i][j][k])
    else:
        return CE*Ex[i][j][k]+CEL*(Hz[i][j][k]-Hz[i][j-1][k])-CEL*(Hy[i][j][k]-Hy[i][j][k-1])
    
def Ey_n(i,j,k):
    if i==0:
        return Ey[i+1][j][k]+mur*(Ey_n(i+1,j,k)-Ey[i][j][k])
    elif i==ny-1:
        return cpEy[i-1][j][k]+mur*(Ey[i-1][j][k]-Ey[i][j][k])
    elif j==0:
        return Ey[i][j+1][k]+mur*(Ey_n(i,j+1,k)-Ey[i][j][k])
    elif j==ny-1:
        return cpEy[i][j-1][k]+mur*(Ey[i][j-1][k]-Ey[i][j][k])
    elif k==0:
        return Ey[i][j][k+1]+mur*(Ey_n(i,j,k+1)-Ey[i][j][k])
    elif k==nz-1:
        return cpEy[i][j][k-1]+mur*(Ey[i][j][k-1]-Ey[i][j][k])
    else:
        return CE*Ey[i][j][k]+CEL*(Hx[i][j][k]-Hx[i][j][k-1])-CEL*(Hz[i][j][k]-Hy[i-1][j][k])

def Ez_n(i,j,k):
    if i==0:
        return Ez[i+1][j][k]+mur*(Ez_n(i+1,j,k)-Ez[i][j][k])
    elif i==ny-1:
        return cpEz[i-1][j][k]+mur*(Ez[i-1][j][k]-Ez[i][j][k])
    elif j==0:
        return Ez[i][j+1][k]+mur*(Ez_n(i,j+1,k)-Ez[i][j][k])
    elif j==ny-1:
        return cpEz[i][j-1][k]+mur*(Ez[i][j-1][k]-Ez[i][j][k])
    elif k==0:
        return Ez[i][j][k+1]+mur*(Ez_n(i,j,k+1)-Ez[i][j][k])
    elif k==nz-1:
        return cpEz[i][j][k-1]+mur*(Ez[i][j][k-1]-Ez[i][j][k])
    else:
        return CE*Ez[i][j][k]+CEL*(Hy[i][j][k]-Hy[i-1][j][k])-CEL*(Hx[i][j][k]-Hx[i][j-1][k])
    
def Hx_n(i,j,k):
    if j==ny-1 and k==nz-1:
        return Hx[i][j][k]-CHL*(0-Ez[i][j][k])+CHL*(0-Ey[i][j][k])
    elif j==ny-1:
        return Hx[i][j][k]-CHL*(0-Ez[i][j][k])+CHL*(Ey[i][j][k+1]-Ey[i][j][k])
    elif k==nz-1:
        return Hx[i][j][k]-CHL*(Ez[i][j+1][k]-Ez[i][j][k])+CHL*(0-Ey[i][j][k])
    else:
        return Hx[i][j][k]-CHL*(Ez[i][j+1][k]-Ez[i][j][k])+CHL*(Ey[i][j][k+1]-Ey[i][j][k])
    
def Hy_n(i,j,k):
    if k==nz-1 and i==nx-1:
        return Hy[i][j][k]-CHL*(0-Ex[i][j][k])+CHL*(0-Ez[i][j][k])
    elif k==nz-1:
        return Hy[i][j][k]-CHL*(0-Ex[i][j][k])+CHL*(Ez[i+1][j][k]-Ez[i][j][k])
    elif i==nx-1:
        return Hy[i][j][k]-CHL*(Ex[i][j][k+1]-Ex[i][j][k])+CHL*(0-Ez[i][j][k])
    else:
        return Hy[i][j][k]-CHL*(Ex[i][j][k+1]-Ex[i][j][k])+CHL*(Ez[i+1][j][k]-Ez[i][j][k])
    
def Hz_n(i,j,k):
    if i==nx-1 and j==ny-1:
        return Hz[i][j][k]-CHL*(0-Ey[i][j][k])+CHL*(0-Ex[i][j][k])
    elif i==nx-1:
        return Hz[i][j][k]-CHL*(0-Ey[i][j][k])+CHL*(Ex[i][j+1][k]-Ex[i][j][k])
    elif j==ny-1:
        return Hz[i][j][k]-CHL*(Ey[i+1][j][k]-Ey[i][j][k])+CHL*(0-Ex[i][j][k])
    else:
        return Hz[i][j][k]-CHL*(Ey[i+1][j][k]-Ey[i][j][k])+CHL*(Ex[i][j+1][k]-Ex[i][j][k])

def gaussian(t):
    to=0.1e-9
    a=(4/to)**2
    return 10*math.pow(math.e,-1*a*(t-to)**2)
cm = plt.cm.jet
cm_list = cm(np.arange(cm.N))
cm_list[:, -1] = 0.3
cm_list[[127,128], -1] = 0
mycmap = ListedColormap(cm_list)    
#Ez[24][0][24]=10
#Ez[24][24][0]=1
#Ez[24][24][24]=1
for n in range(nt+1):
    if n%50==0:
        print(n)
    Time=n*dt
    
#    Ez[int(nx/2)][int(ny/2)][int(nz/2)]+=gaussian(Time)
#    Ez[12][12][12]+=gaussian(Time)
#    Ez[24][24][24]+=gaussian(Time)
#    Ez[15][24][24]+=gaussian(Time)
    Ez[15][15][15]+=gaussian(Time)
    

    cpEz=Ez.copy()
    cpEx=Ex.copy()
    cpEy=Ey.copy()
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                Ex[i][j][k]=Ex_n(i,j,k)
                
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                Ey[i][j][k]=Ey_n(i,j,k)
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                Ez[i][j][k]=Ez_n(i,j,k)
                
    
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                Hx[i][j][k]=Hx_n(i,j,k)
                
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                Hy[i][j][k]=Hy_n(i,j,k)
                
                
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                Hz[i][j][k]=Hz_n(i,j,k)

    if n%30==0:
        print(Time)
        fig=plt.figure()
        ax = fig.add_subplot(111, projection='3d', aspect='equal')
        X,Y,Z = np.meshgrid(range(Ez.shape[0]), range(Ez.shape[1]), range(Ez.shape[2]))
        sc = ax.scatter(X, Y, Z, vmin=-0.5,vmax=0.5,c=Ez, alpha=0.3,marker='o', cmap=mycmap,edgecolor='face')
        fig.patch.set_alpha(0.01)
        fig.colorbar(sc)
        plt.show()
