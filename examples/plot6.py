#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import os

plt.rc('font', size=14)
fig, ax = plt.subplots(figsize=(6,4),dpi=100)

ax.set_xlabel(r'$V_{\rm dc}, \hbar\omega_J/e$')
ax.set_ylabel(r'$I/I_c$')
ax.set_xlim([0,0.8])
ax.set_ylim([0.5,1.8])

filename = 'out'
data = np.loadtxt(filename)
MM = len(data[:,1])
M = int(MM/2)
ax.scatter(data[:M,1], data[:M,0], c='k',s=10)
ax.plot(data[:M,1], data[:M,0], 'k',label='down')
ax.scatter(data[M:,1], data[M:,0], c='r',s=10)
ax.plot(data[M:,1], data[M:,0], 'r',label='up')

ax.legend(loc='best')
plt.grid()
plt.tight_layout()
#plt.savefig("gap_resonance.png",dpi=200) ### Uncomment this to export figure
plt.show()
