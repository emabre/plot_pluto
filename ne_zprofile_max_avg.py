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
import matplotlib.lines as mlines

#plt.close("all")
importlib.reload(prf)
importlib.reload(ut)
importlib.reload(apl)
# plt.ion()

# <codecell>
# Options
sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho4.5e-7-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD-r15'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-7-I90-3.2cmL-1mmD'
measurement = '/home/ema/screenshots/ne_pompili2017_0-30mm_0-12e16.png'
# Capillary half length
l_cap = 1.5e-2
# ---

# The frames of pluto which I want to see (it must be a list of integers)
pluto_nframes = [130]  # [55,80,130]
# Capillary radius
r_cap = 0.5e-3

zlim_plot = 0.5
ne_lim_plot = 1.5e17

z_max = 2.5e-2

#%% Get ne, z-profile
z, ne_max_r, times = apl.ne_avg_over_r(sim, pluto_nframes, 'max')
z, ne_int_r, times = apl.ne_avg_over_r(sim, pluto_nframes, 'integral')
# Reflect data
z = np.concatenate((np.flip(-z, axis=0), z), axis=0)
for ii in range(len(pluto_nframes)):
    ne_int_r[ii] = np.concatenate((np.flip(ne_int_r[ii], axis=0), ne_int_r[ii]), axis=0)
    ne_max_r[ii] = np.concatenate((np.flip(ne_max_r[ii], axis=0), ne_max_r[ii]), axis=0)

# t_I, I = ut.get_currtab(sim)


# <codecell> Plots
fig, ax = plt.subplots(figsize = (5.,3), nrows = len(pluto_nframes))
if len(pluto_nframes) == 1:
    ax = [ax]
ne_meas = mpimg.imread(measurement)
for ii in range(len(pluto_nframes)):
    line1 = ax[ii].plot(z, ne_max_r[ii],
                        label='simulated, '+r'$\mathrm{max}(n_e)$',
                        lw = 1.8,
                        color = 'red')[0]
    line2 = ax[ii].plot(z, ne_int_r[ii],
                        label='simulated, '+r'$\mathrm{mean}_r(n_e)$',
                        lw=1.8,
                        color = 'green')[0]

    ax[ii].imshow(ne_meas,
                  extent=[-1.5, 1.5, 0., 12e16],  # extent=[left,right,bottom,top],
                  aspect='auto',
                  origin='upper'
                  )
    # To make legend
    line3 = mlines.Line2D([], [],
                          color='b',
                          linestyle='-',  # , marker='*', markersize=15,
                          label='measured',
                          lw = 3.)
    ax[ii].legend(handles = [line1, line2, line3])

    ax[ii].set_title('t = {:.0f} ns'.format(times[ii]))
    ax[ii].set_ylabel('Electron density $(\mathrm{cm}^{-3})$')
    ax[ii].set_xlabel('z (cm)')
    # ax[ii].grid()

    # ax[ii].set_xlim([-z_max*1e2, z_max*1e2])
    ax[ii].set_xlim([-l_cap*1e2, l_cap*1e2])
fig.tight_layout()
plt.show()
