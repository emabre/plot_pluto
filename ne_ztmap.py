import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import importlib
import os
import pluto_read_frm as prf
from scipy.constants import mu_0
import utilities as ut
import active_plasma_lens as apl

#plt.close("all")
importlib.reload(prf)
importlib.reload(ut)
importlib.reload(apl)
# plt.ion()

# <codecell>
# Options
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1e5Pa'
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa'
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-rhounif-I90-3.2cm'
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
# ---
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-8-I90-3.2cmL-1mmD'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I90-3.2cmL-1mmD-r60-NTOT8'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-7-I90-3.2cmL-1mmD'
sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD-r15'
# ---

plot_ne_map_each_frame = False

# The frames of pluto which I want to see (it must be a list of integers)
pluto_nframes = list(range(5,290,5))
# z position of z-const lines (in cm)
# Z lines settings, z lines always start from 0
N_z_lines = 30
z_lines_end = 0.5
# Capillary radius
r_cap = 0.5e-3
zlim_plot = 0.5
ne_lim_plot = 1.5e17
average_ne = 'max'  # 'max','integral'
reflect_lowz = True

#%% Get ne, z-profile
z_cc, z, ne_avg_r, t = apl.ne_avg_over_r(sim, pluto_nframes, average_ne,
                                              ret_z_cell_borders=True)
t = np.array(t)
ne_avg_r = np.array(ne_avg_r)
# I average ne over the time cells, so that I get a quantity that is cell centered
# also on the time grid (before it was vertex centered on the time grid and cell centered
# on the z grid)
ne_avg_r_cc = 0.5*(ne_avg_r[1:,:]+ne_avg_r[:-1,:])

if reflect_lowz:
    # Maybe the code would work even without flipping, but I do so, to make if more robust
    z = np.concatenate((np.flip(-z, axis=0), z[1:]))
    ne_avg_r_cc = np.concatenate((np.flip(ne_avg_r_cc, axis=1), ne_avg_r_cc), axis=1)

tt,zz = np.meshgrid(t,z)

#%% Read picture
import matplotlib.image as mpimg
filippi = mpimg.imread('/home/ema/screenshots/filippi.png')
fig, ax = plt.subplots()
ax.imshow(filippi,
          extent=[-1.5, 1.5,1250.,150.],  # extent=[left,right,bottom,top],
          aspect='auto',
          origin='upper'
            )

# <codecell> Plots
fig, ax = plt.subplots()

mp = ax.pcolormesh(zz, tt, ne_avg_r_cc.T,
                   cmap='jet',
                   vmax=1.5e17,
                   vmin=0.0)
ax.invert_yaxis()
ax.set_xlim([-1.5,1.5])
ax.set_ylim([1250., 150.])
fig.colorbar(mp)
