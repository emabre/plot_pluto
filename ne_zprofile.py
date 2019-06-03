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
sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD-r15',
       '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD',
       '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD-r60-NTOT8']
# ---

plot_ne_map_each_frame = False

# The frames of pluto which I want to see (it must be a list of integers)
pluto_nframes = list(range(5,120,30))
# z position of z-const lines (in cm)
# Z lines settings, z lines always start from 0
N_z_lines = 30
z_lines_end = 0.5
# Capillary radius
r_cap = 0.5e-3

reflect_lowz = True
zlim_plot = 0.5
ne_lim_plot = 1.5e17

average_ne = 'max'  # 'max','integral'

#%% Get ne, z-profile
ne_avg_r = [[],[],[]]
for ss in range(len(sim)):
    z, ne_avg_r[ss], times = apl.ne_avg_over_r(sim[ss], pluto_nframes, average_ne)
# t_I, I = ut.get_currtab(sim)

# <codecell> Plots
fig, ax = plt.subplots(nrows=len(pluto_nframes), sharex=True, sharey=True)

for ff in range(len(pluto_nframes)):
    for ss in range(len(sim)):
        label = os.path.basename(sim[ss])
        if label=='':  # This is in case the sim path ends with '/'
            label = os.path.basename(sim[ss][:-1])
        ax[ff].plot(z, ne_avg_r[ss][ff], label=label)
    ax[ff].set_title('t = {} ns'.format(times[ff]))
    ax[ff].set_ylabel('ne (cm^-3)')
ax[-1].set_xlabel('z (cm)')
ax[0].legend()
plt.tight_layout()
plt.show()
