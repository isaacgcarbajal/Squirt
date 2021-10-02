import numpy as np
import matplotlib.pyplot as plt
import sys

fname = sys.argv[1]

neq = 4
nx = 256
ny = 256

x = np.linspace(0,2,nx+1)
y = np.linspace(0,2,ny+1)

data = np.fromfile(fname, dtype="d", count=neq*nx*ny).reshape((neq,nx,ny))

rho = data[0,:,:]
u = data[1,:,:]
v = data[2,:,:]
P = data[3,:,:]

fig, ax  = plt.subplots(2,2, figsize=(9, 9))

cm = ax[0,0].imshow(rho.T, origin='lower'); fig.colorbar(cm,ax=ax[0,0])
cm = ax[0,1].imshow(u.T, origin='lower');   fig.colorbar(cm, ax=ax[0,1])
cm = ax[1,0].imshow(v.T, origin='lower');   fig.colorbar(cm, ax=ax[1,0])
cm = ax[1,1].imshow(P.T, origin='lower');   fig.colorbar(cm, ax=ax[1,1])

figNum=int(fname.split('_')[1].split('.')[0])
figTitle=r"$t = %d~[\rm{dt^{-1}}]$"%(figNum)

dy = 0.2
ax[0,0].set(
        xlabel=r"$x$",
        ylabel=r"$y$",
        title=r"$\rho(x)$"
        )

ax[0,1].set(
        xlabel=r"$x$",
        ylabel=r"$y$",
        title=r"$u(x)$"
        )

ax[1,0].set(
        xlabel=r"$x$",
        ylabel=r"$y$",
        title=r"$v(x)$",
        )

ax[1,1].set(
        xlabel=r"$x$",
        ylabel=r"$y$",
        title=r"$P(x)$"
        )


fig.text(0.45,0.95,figTitle)



fig.savefig(fname.split('.')[0]+'.png')
