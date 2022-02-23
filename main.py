import csv
import numpy as np
import matplotlib.pyplot as plt

# function imports
from vs_polar import vs_Efield_polar2D
from vs_polar import vs_Efield_polar2D_M2
from vs_testdata import vs_Efield_testdata

# constants
RE      = 6378.1370  # equatorial radius of Earth; [km]
gamma   = 2          # shielding factor
kp      = 1          # kp index

rmax = 3*RE
dn = 300
'''
Calculations
'''
# METHOD 1: using numpy gradient to calculate E
r, phi, edat1 = vs_Efield_polar2D(rmax, dn, gamma, kp)

# METHOD 2: using derived components of E (Er, Ep)
r, phi, edat2 = vs_Efield_polar2D_M2(rmax, dn, gamma, kp)

'''
Plots
'''
pcolor = 'gnuplot2'

fig, axes = plt.subplots(nrows=1, ncols=2, squeeze=False, subplot_kw=dict(projection='polar'))
fig.suptitle("Volland-Stern Electric Field $E_c$ ($K_p =$" + str(kp) + ")")

# M1
ax1 = axes[0,0]
ax1.set_title('Method 1 (np gradient)')

fig1 = ax1.pcolormesh(phi,r,edat1, cmap=pcolor, shading='nearest')
fig.colorbar(fig1, ax=ax1)
ax1.grid(color='w', linestyle=':', linewidth=0.75)

# M2
ax2 = axes[0,1]
ax2.set_title('Method 2 (derived E)')
fig2 = ax2.pcolormesh(phi,r,edat2, cmap=pcolor, shading='nearest')
fig.colorbar(fig2, ax=ax2)
ax2.grid(color='w', linestyle=':', linewidth=0.75)


plt.show()