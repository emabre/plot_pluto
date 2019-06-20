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
sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD-NEWGRID',
       '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD-NEWGRID'
       ]
names = ('(a)', '(b)')
# ---

plot_ne_map_each_frame = False

# The frames of pluto which I want to see (it must be a list of integers)
pluto_nframes = [10, 24, 40, 80] # list(range(5,120,30))
# z position of z-const lines (in cm)
# Z lines settings, z lines always start from 0
N_z_lines = 30
z_lines_end = 0.5
# Capillary radius
r_cap = 0.5e-3
l_cap = 0.5e-2

reflect_lowz = True
zlim_plot = 0.5
ne_lim_plot = 1.5e17

average_ne = 'axis'  # 'axis', 'max','integral'

#%% Get ne, z-profile
ne_avg_r = [[],[],[]]
for ss in range(len(sim)):
    z, ne_avg_r[ss], times = apl.ne_avg_over_r(sim[ss], pluto_nframes, average_ne)
# t_I, I = ut.get_currtab(sim)

# <codecell> Plots
if len(pluto_nframes) != 4:
    fig, ax = plt.subplots(nrows=len(pluto_nframes), sharex=True, sharey=True)

    for ff in range(len(pluto_nframes)):
        for ss in range(len(sim)):
            ax[ff].plot(z, ne_avg_r[ss][ff], label=names[ss])
        ax[ff].set_title('t = {} ns'.format(times[ff]))
        ax[ff].set_ylabel('ne (cm^-3)')
    ax[-1].set_xlabel('z (cm)')
    ax[0].legend()
    plt.tight_layout()
    plt.show()
else:
    fig, ax = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True, figsize=(5.5, 4.5))
    idx = [(0,0), (0,1), (1,0), (1,1)]
    for ff in range(len(pluto_nframes)):
        for ss in range(len(sim)):
            ax[idx[ff]].semilogy(z, ne_avg_r[ss][ff], label=names[ss])
        ax[idx[ff]].set_title('t = {:.0f}ns'.format(times[ff]))
        # ax[idx[ff]].set_ylabel('Electron density $(\mathrm{cm}^{-3})$')
        ax[idx[ff]].axvspan(0., l_cap*1e2, facecolor='gray', alpha=0.4)
        ax[idx[ff]].grid()
    ax[0,0].set_xlim([0., 2.])
    ax[0,0].set_ylim([1e10, 1e18])
    ax[-1,0].set_xlabel('z (cm)')
    ax[-1,1].set_xlabel('z (cm)')
    for ii in (0,1):
        ax[ii,0].set_ylabel('Electron density $(\mathrm{cm}^{-3})$')
    ax[0,0].legend()
    plt.tight_layout()
    plt.show()
