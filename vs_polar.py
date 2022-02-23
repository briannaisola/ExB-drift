import numpy as np
import matplotlib.pyplot as plt

# METHOD 1: using numpy gradient to calculate E
def vs_potential(r_geo, phi, gamma=2, kp=1):
    '''
    DESC.:  compute vs potential in 2D spherical coordinates
    INPUT:  r_geo       ; radial geocentric distance
            phi       ; azimuthal angle, MLT with respect to noon
            gamma*      ; shielding parameter   ; DEFAULT = 2
            kp*         ; kp index  ; DEFAULT = 1
    '''

    RE = 6378.1370  # equatorial radius of Earth; [km]
    E0 = 0.045 / ((1 - (0.159 * kp) + (0.0093 * kp ** 2)) ** 3 * (RE ** 2))

    U = -E0*(r_geo**gamma)*np.sin(phi)
    return U

def vs_Efield_polar2D(rmax, dn, gamma=2, kp=1):
    '''
    DESC.:  compute VS electric field in 2D spherical coordinates
    INPUT:  rmax        ; max radial geocentric distance
            dn          ; number of grid steps
            gamma*      ; shielding parameter   ; DEFAULT = 2
            kp*         ; kp index  ; DEFAULT = 1
    OUTPUT: dr x dr array of electric field values
    '''

    angle = 2.0 * np.pi  # compute around 360-deg

    rdat = np.linspace(0, rmax, dn)
    tdat = np.linspace(0.0, angle, dn)

    r, phi = np.meshgrid(rdat, tdat)

    # array to hold potential field
    umap = np.zeros(r.shape)

    # calculate potential map
    for i in range(0, dn):
        for j in range(0, dn):
            rval = r[i, j]
            tval = phi[i, j]

            # asign values
            umap[i, j] = vs_potential(rval, tval, gamma, kp)

    # calculate electric field using numpy gradient
    egrad = np.gradient(umap)
    emap = np.sqrt(egrad[0] ** 2 + egrad[1] ** 2)

    return r, phi, emap

# METHOD 2: using derived components of E (Er, Ep)
def vs_Efield_polar2D_M2(rmax, dn, gamma=2, kp=1):
    '''
    DESC.:  compute VS electric field in 2D spherical coordinates
    INPUT:  rmax        ; max radial geocentric distance
            dn          ; number of grid steps
            gamma*      ; shielding parameter   ; DEFAULT = 2
            kp*         ; kp index  ; DEFAULT = 1
    OUTPUT: dr x dr array of electric field values
    '''

    angle = 2.0 * np.pi  # compute around 360-deg

    rdat = np.linspace(0, rmax, dn)
    tdat = np.linspace(0, angle, dn)

    r, phi = np.meshgrid(rdat, tdat)

    # array to hold electric field
    emap = np.zeros(r.shape)

    # calculate potential map
    RE = 6378.1370  # equatorial radius of Earth; [km]
    E0 = 0.045 / ((1 - (0.159 * kp) + (0.0093 * kp ** 2)) ** 3 * (RE ** 2))

    for i in range(0, dn):
        for j in range(0, dn):
            rval = r[i, j]
            tval = phi[i, j]

            # compute F-field components
            Er = -E0 * gamma * (rval ** (gamma - 1)) * np.sin(tval)
            Ep = -E0 * (rval ** (gamma - 1)) * np.cos(tval)

            # test
            

            # write magnitude to emap
            emap[i,j] = np.sqrt(Er ** 2 + Ep ** 2)

    return r, phi, emap