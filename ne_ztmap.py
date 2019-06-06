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
# measurement = 'file:///home/ema/Dottorato/dati_sperimentali_e_calcoli/misure_ne_da_confrontare/misure_capillare_1cmL-1mmD/stessa_scala_temporale/ne_I245_t50-1450ns_z0-10mm_cmapjet0-15e16_CUT.png'
measurement = '/home/ema/Dottorato/dati_sperimentali_e_calcoli/misure_ne_da_confrontare/misure_capillare_1cmL-1mmD/stessa_scala_temporale/ne_I245_t50-1450ns_z0-10mm_cmapjet0-15e16_CUT.png'

# ---
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-8-I90-3.2cmL-1mmD'
sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I90-3.2cmL-1mmD-r60-NTOT8'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho4.5e-7-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-7-I90-3.2cmL-1mmD'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD-r15'
measurement = '/home/ema/Dottorato/dati_sperimentali_e_calcoli/misure_ne_da_confrontare/misure_capillare_3cmL-1mmD/stessa_scala_temporale/densità_elettronica/ne_I90_t150-1250ns_z0-30mm_cmapjet0-15e16_CUT.png'
# ---

# The frames of pluto which I want to see (it must be a list of integers)
pluto_nframes = list(range(0,221,10))
# Capillary radius
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
gs = gridspec.GridSpec(1, 3,
                       width_ratios=[28, 28, 2],
                       # height_ratios=[1, 1]
                       )
fig = plt.figure(figsize=(4,3))
ax_sim  = plt.subplot(gs[0,0])
ax_meas = plt.subplot(gs[0,1])
ax_cb = plt.subplot(gs[0,2])

filippi = mpimg.imread(measurement)
# fig, ax = plt.subplots(ncols=2)
ax_meas.imshow(filippi,
             extent=[-1.5, 1.5,1250.,150.],  # extent=[left,right,bottom,top],
             aspect='auto',
             origin='upper'
             )
ax_meas.set_yticks([],[])
mp = ax_sim.pcolormesh(zz, tt, ne_avg_r_cc.T,
                   cmap='jet',
                   vmax=1.5e17,
                   vmin=0.0)

ax_sim.invert_yaxis()
ax_sim.set_xlim([-1.5,1.5])
ax_sim.set_ylim([1250., 150.])
xticks = [-1.5,0.,1.5]
ax_meas.set_xticks(xticks)
ax_sim.set_xticks(xticks)

ax_sim.set_ylabel('Time (s)')
ax_sim.set_xlabel('Longitudinal position (mm)')
ax_sim.set_xlabel('Longitudinal position (mm)')
ax_sim.set_title('Simulation')
ax_meas.set_title('Measurement')

fig.colorbar(mp, cax = ax_cb)
fig.tight_layout()
