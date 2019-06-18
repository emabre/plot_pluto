import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import importlib
import os
import pluto_read_frm as prf
from scipy.constants import mu_0
import utilities as ut
import active_plasma_lens as apl
from scipy.interpolate import griddata
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
from copy import copy

#plt.close("all")
importlib.reload(prf)
importlib.reload(ut)
importlib.reload(apl)

# <codecell>
# Settings
# Available: 'compare-1mmD-1.2mmD', 'compare-pI500-mI500'
setting_case = 'compare-1mmD-1.2mmD'

# ----
sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD'
pluto_nframe = 10  # 10
# Capillary length, half of the real one (including electrodes)
l_cap = 0.6e-2  # m
r_cap = 0.5e-3  # m
dz_cap = 1.e-3  # m
rmax = 1e-3  # m
zmax = 1.2e-2  # m

#%% Load the data
pluto_dir = os.path.join(os.path.expandvars(sim), 'out')
q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir, nframe=pluto_nframe)
T = q["T"]
B = q["bx3"]*1e-4  # Convert to Tesla
vr = q["vx1"]*1e-2  # Converto to m/s
vz = q["vx2"]*1e-2  # Converto to m/s
rho = q["rho"]*1e-3  # Convert to g/cm3
# Convert r and z to m
r /= 1e5
z /= 1e5

#%%
# Compute the cell centers
r_cc = 0.5*(r[1:]+r[:-1])
z_cc = 0.5*(z[1:]+z[:-1])

#%% Plotting
gs = gridspec.GridSpec(6, 2,
                       width_ratios=[30, 2],
                       height_ratios=6*[1],
                       hspace = 0
                       )

fig = plt.figure(figsize=(5.5,4.5))
ax = [[],[]]; ax_cb = [[],[]]
ax[0] = plt.subplot(gs[:3,0])
ax[1] = plt.subplot(gs[3:,0])
ax_cb[0] = plt.subplot(gs[:2,1])
ax_cb[1] = plt.subplot(gs[4:,1])

mp_T = ax[0].contourf(z_cc*1e2, r_cc*1e6,
                    rho.T,
                     # levels = np.linspace(Tmin,Tmax,101),
                     # levels = 100,
                     # levels=np.linspace(Tmin, Tmax, 150),
                     # locator=ticker.LogLocator(),
                     # norm = colors.LogNorm(vmin=vmin, vmax=vmax),
                     # vmin = Tmin,
                     # vmax = Tmin,
                     # cmap = 'hot',
                     cmap = 'pink',
                     )

# Tmin = 2.5e3
# Tmax = 80e3
# for ss in (0,1):
#     T[ss][T[ss]<2.5e3] = 2.5e3
# for ss in (0,1):
mp_rho = ax[1].contourf(z_cc*1e2, r_cc*1e6,
                     T.T,
                     # levels = np.linspace(Tmin,Tmax,101),
                     # levels = 100,
                     # levels=np.linspace(Tmin, Tmax, 150),
                     # locator=ticker.LogLocator(),
                     # norm = colors.LogNorm(vmin=vmin, vmax=vmax),
                     # vmin = Tmin,
                     # vmax = Tmin,
                     # cmap = 'hot',
                     cmap = 'gist_heat',
                     )
    # # Plot iso-I. Note that in cylindrical coordinates, their density does not represent current density (even though correlated),
    # # because the current density is not divergent free in cylindrical symmetry with the cartesia definition of divergence ((Dz,Dr)).
    # mpj = ax[ss].contour(z_cc[ss]*1e2, r_cc[ss]*1e6,
    #                      (B[ss]*r_cc[ss]*2*np.pi/mu_0).T,  # from Ampere's law, I=B*r*2*pi/mu_0
    #                      levels = 7,
    #                      colors='cyan',
    #                      linewidths = 1.,
    #                     )

ax[0].set_ylim(0.01, rmax*1e6)
ax[1].set_ylim(rmax*1e6, 0.01)
for ii in (0,1):
    ax[ii].set_xlim(0.01, zmax*1e2)
ax[0].set_xticklabels([])

ax[1].set_xlabel('z (cm)')
ax[0].set_ylabel('r (μm)')
ax[1].set_ylabel('r (μm)')

cbar_T = fig.colorbar(mp_T, cax = ax_cb[0],
                    label='Temperature (K)',
                    # ticks = [2500, 2e4, 4e4, 6e4, 8e4],
                    # ticks = [2500, 1e4, 2e4, 3e4, 4e4],
                    )
cbar_rho = fig.colorbar(mp_rho, cax = ax_cb[1],
                    label='Mass density (g/cm3)',
                    # ticks = [2500, 2e4, 4e4, 6e4, 8e4],
                    # ticks = [2500, 1e4, 2e4, 3e4, 4e4],
                    )
# cbar.set_ticklabels(['$<10^8$', '$10^{10}$', '$10^{12}$', '$10^{14}$', '$10^{16}$', '$10^{18}$'])

# ----
# Draw electrodes and capillary walls
# Locations of lower left corners of electrodes z(m),r(m) (one for each sim)
electrode_zr_ll = ((l_cap-dz_cap)*1e2, r_cap*1e6)
wall_zr_ll = (0, r_cap*1e6)

electrodes = [Rectangle((electrode_zr_ll[0], electrode_zr_ll[1]), dz_cap*1e2, r.max()*1e6,
                        fill=True, facecolor='#a3552c', edgecolor='k', zorder=10 ) for ss in (0,1)]
walls = [Rectangle((wall_zr_ll[0], wall_zr_ll[1]), (l_cap-dz_cap)*1e2, r.max()*1e6,
                        fill=True, facecolor='gray', edgecolor='k', zorder=10 ) for ss in (0,1)]
for ss in (0,1):
    ax[ss].add_patch(electrodes[ss])
    ax[ss].add_patch(walls[ss])
# el_coll = PatchCollection([el], facecolor='r',
#                           alpha = 0.99,
#                           edgecolor='g')
# ax[0].add_collection(el_coll)

# ----
ax[0].set_title('t={:g}ns'.format(t))

fig.tight_layout()
plt.show()
