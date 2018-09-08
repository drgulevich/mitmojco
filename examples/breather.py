#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

data = np.loadtxt("breather.dat")
frames=data.shape[0]
N=data.shape[1]
L0=40. # normalized length of the Josephson junction
x=np.linspace(-L0/2.,L0/2.,N)
print('L0:',L0)

plt.rc('font', family='serif', size='14')
fig, ax = plt.subplots(figsize=(8,6))
ax.set_xlim([-L0/2.,L0/2.])
ax.set_ylim([-7,7])

ax.set_title("Sine-Gordon breather")
ax.set_xlabel(r'$x,\lambda_J$',fontsize=18)
ax.set_ylabel(r'$\varphi(x,t)$',fontsize=18)
line, = ax.plot([], [], lw=2)
time_template = 'Frame number: %d'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(i):
    y=data[i,:]
    line.set_data(x, y)
    time_text.set_text(time_template%(i))
    return line, time_text

anim = animation.FuncAnimation(fig, animate, init_func=init, frames=frames, interval=100, blit=True, repeat=False)

#ax.grid()
plt.show()
