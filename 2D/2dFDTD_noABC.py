import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import os
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
import datetime

# 時刻によりファイルを作成

simulation_time = datetime.datetime.today()
folder_name = "2DFDTD_noABC_{0:%Y%m%d-%H%M%S}".format(simulation_time)
os.mkdir(folder_name)


def gausiannpalse(t):
    to = 0.5e-9
    a = (4 / to) ** 2
    return math.pow(math.e, -1 * a * (t - to) ** 2)


fig = plt.figure()
# 初期化、今回は異物を空気とする
step_i = 161
step_j = 201
timestep = 2000
anomary_step_i = [64, 96]
anomary_step_j = [80, 120]

H1 = np.zeros((step_i, step_j))
H2 = np.zeros((step_i, step_j))
E3 = np.zeros((step_i, step_j))
ID = np.zeros((step_i, step_j))
params_ = np.zeros((step_i, step_j, 3))
id = np.zeros((2, 3))
params = np.zeros((2, 3))
dx = 2.0e-3  # f=250,dx=λ*0.25
dt = 3.0e-12
"""
骨材とコンクリートの設定
"""
for i in range(2):
    if i == 0:  # コンクリート id:0
        ε = 66.41 * (1e-12)
        μ = 1.257 * (1e-6)
        σ = 0.001
    if i == 1:  # 骨材 id:1
        ε = 8.855 * (1e-12)
        μ = 1.257 * (1e-6)
        σ = 0.000000000001
    params[i][0] = dt / (μ * dx)
    params[i][1] = (1 - σ * dt / 2 / ε) / (1 + σ * dt / 2 / ε)
    params[i][2] = (dt / (dx * ε)) / (1 + σ * dt / 2 / ε)

# 一種類、一個のみ
params_[:] = np.copy(params[0])
for i in range(anomary_step_i[0], anomary_step_i[1]):
    for j in range(anomary_step_j[0], anomary_step_j[1]):
        params_[i][j] = np.copy(params[1])

ims = []
# 電磁場の計算

for n in range(timestep + 1):
    print(n)
    T = n * dt
    E3[80][100] += gausiannpalse(T)

    E3[0,0] = params_[i][0,0] * E3[0,0] + params_[0,0,2] * (H2[0,0] - 0 - H1[0,0] + 0)
    E3[0,1:] = params_[0,1:,1] * E3[0,1:] + params_[0,1:,2] * (H2[0,1:] - 0 - H1[0,1:] + H1[0,0:-1])
    E3[1:,0] = params_[1:,0,1] * E3[1:,0] + params_[1:,0,2] * (H2[1:,0] - H2[0:-1,0] - H1[1:,0] + 0)
    E3[1:,1:] = params_[1:,1:,1] * E3[1:,1:] + params_[1:,1:,2] * (H2[1:,1:] - H2[0:-1,1:] - H1[1:,1:] + H1[1:,0:-1])

    H2[-1,0:-1] + params_[-1,0:-1,0] * (0 - E3[-1,0:-1])
    H2[0:-1,0:-1] = H2[0:-1,0:-1] + params_[0:-1,0:-1,0] * (E3[1:,0:-1] - E3[0:-1,0:-1])

    H1[0:-1,-1] = H1[0:-1,-1] - params_[0:-1,-1,0] * (0 - E3[0:-1,-1])
    H1[0:-1,0:-1] = H1[0:-1,0:-1] - params_[0:-1,0:-1,0] * (E3[0:-1,1:] - E3[0:-1,0:-1])

    if n % 20 == 0:
        """
        #静止画表示
        ax = fig.add_subplot(111)
        cax = ax.imshow(E3, interpolation='none',vmin=-0.5,vmax=0.5,animated=True)
        fig.colorbar(cax)
        """
        im = plt.imshow(E3, vmin=-0.05, vmax=0.05, animated=True)
        ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=50,
                                repeat_delay=1000)

plt.show()
