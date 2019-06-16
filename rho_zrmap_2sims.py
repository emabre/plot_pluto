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
if setting_case == 'compare-1mmD-1.2mmD':
    sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD',
           '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD'
           ]
    # Capillary length, half of the real one (including electrodes)
    l_cap = 0.6e-2  # m
    r_cap = (0.5e-3, 0.6e-3)  # m
    dz_cap = 1.e-3  # m
    rmax = 2e-3  # m
    zmax = 1.4e-2  # m
    pluto_nframe = 0
    # Scale choice: 'lin' or 'log'
    scale = 'lin'

elif setting_case == 'compare-pI500-mI500':
    sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-mI550-3.2cmL-1mmD-r60-NTOT20-diffRecPeriod10-NB20/',
           '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I550-3.2cmL-1mmD-r60-NTOT32-diffRecPeriod8/'
           ]
    pluto_nframe = 60  # 10
    # Capillary length, half of the real one (including electrodes)
    l_cap = 1.6e-2  # m
    r_cap = (0.5e-3, 0.5e-3)  # m
    dz_cap = 1.e-3  # m
    rmax = 2.e-3  # m
    zmax = 2.8e-2  # m
    # Scale choice: 'lin' or 'log'
    scale = 'log'

#%% Load the data
if len(sim)!=2:
    raise ValueError('sim must be a list of two(2) simulation paths')
rho = [[], []]
r = [[],[]]
z = [[],[]]
t = [[],[]]
for ss in range(len(sim)):
    pluto_dir = os.path.join(os.path.expandvars(sim[ss]), 'out')
    q, r[ss], z[ss], theta, t[ss], n = prf.pluto_read_vtk_frame(pluto_dir,
                                                                # time=125.0,
                                                                nframe=pluto_nframe)
    rho[ss] = q["rho"]
    # Convert r and z to m
    r[ss] /= 1e5
    z[ss] /= 1e5

#%%
# Compute the cell centers
r_cc = [0.5*(r[ss][1:]+r[ss][:-1]) for ss in (0,1)]
z_cc = [0.5*(z[ss][1:]+z[ss][:-1]) for ss in (0,1)]

#%% Plotting
gs = gridspec.GridSpec(2, 2,
                       width_ratios=[30, 2],
                       height_ratios=[1, 1]
                       )

fig = plt.figure(figsize=(5.5,4.5))
ax = [[],[]]
ax[0] = plt.subplot(gs[0,0])
ax[1] = plt.subplot(gs[1,0])
ax_cb = plt.subplot(gs[:,1])
if scale=='log':
    rho_logmin = -9
    rho_logmax = -6
    for ss in (0,1):
        rho[ss][rho[ss] < 10**rho_logmin] = 10**rho_logmin
        mp = ax[ss].contourf(z_cc[ss]*1e2, r_cc[ss]*1e6,
                             # rho[ss].T,
                             np.log10(rho[ss].T),
                             levels=np.linspace(rho_logmin, rho_logmax, 150),
                             # locator=ticker.LogLocator(),
                             # norm = colors.LogNorm(vmin=vmin, vmax=vmax),
                             cmap = 'pink')
elif scale=='lin':
    rho_min = 2.5e-10
    rho_max = 2.5e-7
    for ss in (0,1):
        rho[ss][rho[ss] < rho_min] = rho_min
        mp = ax[ss].contourf(z_cc[ss]*1e2, r_cc[ss]*1e6,
                             # rho[ss].T,
                             rho[ss].T,
                             levels = np.linspace(0, rho_max, 151),
                             # vmin = rho_min,
                             # vmax = rho_max,
                             # norm = colors.LogNorm(vmin=vmin, vmax=vmax),
                             cmap = 'pink')
else:
    raise ValueError('Wrong setting for scale')

if scale=='log':
    ticks = list(range(rho_logmin, rho_logmax+1, 1))
    cbar = fig.colorbar(mp, cax = ax_cb,
                 label='Mass density $(\mathrm{g} / \mathrm{cm}^{-3})$',
                 ticks = ticks)
    cbar.set_ticklabels(['$<10^{{{}}}$'.format(ticks[0])] + ['$10^{{{}}}$'.format(ticks[ii]) for ii in range(1,len(ticks))])
elif scale=='lin':
    ticks = [2.5e-10, 0.5e-7, 1.e-7, 1.5e-7, 2.e-7, 2.5e-7]
    cbar = fig.colorbar(mp, cax = ax_cb,
                        label='Mass density $(\mathrm{g} / \mathrm{cm}^{-3})$',
                        ticks = ticks,
                        )
    # cbar.set_ticklabels(['$<10^{{{}}}$'.format(ticks[0])] + ['$10^{{{}}}$'.format(ticks[ii]) for ii in range(1,len(ticks))])

#---
# Labels and lims
for ss in (0,1):
    ax[ss].set_ylim(0.01, rmax*1e6)
    ax[ss].set_xlim(0., zmax*1e2)
ax[1].set_xlabel('z (cm)')
ax[0].set_ylabel('r (μm)')
ax[1].set_ylabel('r (μm)')

# ----
# Draw electrodes and capillary walls
# Locations of lower left corners of electrodes z(m),r(m) (one for each sim)
electrode_zr_ll = (((l_cap-dz_cap)*1e2, r_cap[0]*1e6),  # sim 0
                   ((l_cap-dz_cap)*1e2, r_cap[1]*1e6))  # sim 1
wall_zr_ll = ((0, r_cap[0]*1e6),  # sim 0
              (0, r_cap[1]*1e6))  # sim 1

electrodes = [Rectangle((electrode_zr_ll[ss][0], electrode_zr_ll[ss][1]), dz_cap*1e2, r[ss].max()*1e6,
                        fill=True, facecolor='#a3552c', edgecolor='k' ) for ss in (0,1)]
walls = [Rectangle((wall_zr_ll[ss][0], wall_zr_ll[ss][1]), (l_cap-dz_cap)*1e2, r[ss].max()*1e6,
                        fill=True, facecolor='gray', edgecolor='k' ) for ss in (0,1)]
for ss in (0,1):
    ax[ss].add_patch(electrodes[ss])
    ax[ss].add_patch(walls[ss])
# el_coll = PatchCollection([el], facecolor='r',
#                           alpha = 0.99,
#                           edgecolor='g')
# ax[0].add_collection(el_coll)

# ----

if setting_case == 'compare-1mmD-1.2mmD':
    case = ('(a)', '(b)')
    for ss in (0,1):
        ax[ss].set_title('{:s}, t={:g}ns, R={:g}μm'.format(case[ss], t[0], r_cap[ss]*1e6))
else:
    for ss in (0,1):
        ax[ss].set_title('time: {:g}ns'.format(t[0]))
fig.tight_layout()
plt.show()
# rr = [[],[]]
# zz = [[],[]]
# for ss in (0,1):
#     rr[ss], zz[ss] = np.meshgrid(r[ss], z[ss])
#
# fig, ax = plt.subplots()
# ss = 0
# ax.pcolormesh(rr[0], zz[0], ne_cc[0])
# ax.pcolormesh(rr[0], zz[0], ne[0])
