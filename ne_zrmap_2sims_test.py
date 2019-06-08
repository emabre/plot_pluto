import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import importlib
import os
import pluto_read_frm as prf
from scipy.constants import mu_0
import utilities as ut
import active_plasma_lens as apl
import matplotlib.image as mpimg
from scipy.interpolate import griddata

#plt.close("all")
importlib.reload(prf)
importlib.reload(ut)
importlib.reload(apl)

# <codecell>
# Settings

# ----
sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD',
       '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD'
       ]
pluto_nframe = 50
l_cap = 1e-2 # m
r_cap = 0.5e-2

#%% Load the data
if len(sim)!=2:
    raise ValueError('sim must be a list of two(2) simulation paths')
ne = [[], []]
r = [[],[]]
z = [[],[]]
t = [[],[]]
for ss in range(len(sim)):
    pluto_dir = os.path.join(os.path.expandvars(sim[ss]), 'out')
    q, r[ss], z[ss], theta, t[ss], n = prf.pluto_read_vtk_frame(pluto_dir,
                                                                # time=125.0,
                                                                nframe=pluto_nframe)
    ne[ss] = q["ne"]
    # Convert r and z to m
    r[ss] /= 1e5
    z[ss] /= 1e5

#%%
# Interpolate with nearest the border points

r_cc = [0.5*(r[ss][1:]+r[ss][:-1]) for ss in (0,1)]
z_cc = [0.5*(z[ss][1:]+z[ss][:-1]) for ss in (0,1)]

# fig, ax = plt.subplots(nrows=2)

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

ss = 0
ax[0].contourf(z_cc[0]*1e2, r_cc[0]*1e6, ne[0].T,
            levels=100,
            cmap = 'inferno')
mp = ax[1].contourf(z_cc[1]*1e2, r_cc[1]*1e6, ne[1].T,
                    levels=100,
                    cmap = 'inferno')

# ax[0].set_xticks([])

ax[1].set_ylim(0., 2e3)
ax[0].set_ylim(0.,2e3)
for ss in (0,1):
    ax[ss].set_xlim(0.01,1.6)
ax[1].set_xlabel('z (cm)')
ax[0].set_ylabel('r (μm)')
ax[1].set_ylabel('r (μm)')

fig.colorbar(mp, cax = ax_cb, label='Electron density $(\mathrm{cm}^{-3})$')
fig.tight_layout()
# rr = [[],[]]
# zz = [[],[]]
# for ss in (0,1):
#     rr[ss], zz[ss] = np.meshgrid(r[ss], z[ss])
#
# fig, ax = plt.subplots()
# ss = 0
# ax.pcolormesh(rr[0], zz[0], ne_cc[0])
# ax.pcolormesh(rr[0], zz[0], ne[0])
