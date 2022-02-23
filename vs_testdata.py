import csv
import numpy as np
import matplotlib.pyplot as plt
from vs_polar import vs_potential

# constants
RE      = 6378.1370  # equatorial radius of Earth; [km]
gamma   = 2          # shielding factor
kp      = 1          # kp index


# USING TEST DATA?
# open and read file
file = open('sampledata_9_10_15_to_9_16_2016.csv')
fdat = csv.DictReader(file)

Ldat = []
MLTdat = []

# assign cols to arrays
for col in fdat:
    Ldat.append(float(col['L']))
    MLTdat.append(float(col['MLT']))

# From calculate_potential.py: nL is the number of L values in E, not Φ -- not sure if I did this right
nL = int((max(Ldat)-min(Ldat) + 1))
nMLT = 24

# set and reshape data to arrays
phi = np.array(MLTdat)
phidat = phi.reshape(nL, nMLT)
r = np.array(Ldat)
rdat = r.reshape(nL, nMLT)


def vs_Efield_testdata(r, phi, gamma=2, kp=1):
    size = r.shape
    #print(size)

    # Method 1: calculate potential map and take the gradient
    umap = np.zeros(r.shape)

    size = r.shape

    for i in range(0, size[0]):
        for j in range(0, size[1]):
            rval = r[i, j]
            tval = phi[i, j]

            # asign values
            umap[i, j] = vs_potential(rval, tval, gamma, kp)

    # calculate electric field (take gradient)
    egrad = np.gradient(umap)
    emap_M1 = np.sqrt(egrad[0] ** 2 + egrad[1] ** 2)

    # Method 2: Use E-components
    emap_M2 = np.zeros(r.shape)
    E0 = 0.045 / ((1 - (0.159 * kp) + (0.0093 * kp ** 2)) ** 3 * (RE ** 2))

    for i in range(0, size[0]):
        for j in range(0, size[1]):
            rval = r[i, j]
            tval = phi[i, j]

            Er = E0 * gamma * (rval ** (gamma - 1)) * np.sin(tval)
            Ep = E0 * (rval ** (gamma - 1)) * np.cos(tval)
            emap_M2[i, j] = np.sqrt(Er ** 2 + Ep ** 2)

    return emap_M1, emap_M2


em1, em2 = vs_Efield_testdata(rdat, phidat, gamma=2, kp=1)

fig, axes = plt.subplots(nrows=1, ncols=2, squeeze=False, subplot_kw=dict(projection='polar'))
pcolor = 'gnuplot2'

# M1
ax1 = axes[0,0]

fig1 = ax1.pcolormesh(phidat,rdat,em1, cmap=pcolor, shading='nearest')
fig.colorbar(fig1, ax=ax1)
ax1.grid(color='w', linestyle=':', linewidth=0.75)

# M2
ax2 = axes[0,1]
fig2 = ax2.pcolormesh(phidat,rdat,em2, cmap=pcolor, shading='nearest')
fig.colorbar(fig2, ax=ax2)
ax2.grid(color='w', linestyle=':', linewidth=0.75)





plt.show()