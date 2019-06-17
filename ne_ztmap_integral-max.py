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
# Settings

# Allowed choices are '245A-1cm' and '90A-3cm'
measure_choice = '90A-3cm'  # '90A-3cm' or '245A-1cm'

if measure_choice=='245A-1cm':
    sim = {'mean': '/home/ema/simulazioni/sims_pluto/perTesi/rho6e-7-I245-1.2cmL-1mmD-NEWGRID',
           'max' : '/home/ema/simulazioni/sims_pluto/perTesi/rho2.9e-7-I245-1.2cmL-1mmD-NEWGRID'}
    measurement = '/home/ema/Dottorato/dati_sperimentali_e_calcoli/misure_ne_da_confrontare/misure_capillare_1cmL-1mmD/stessa_scala_temporale/ne_I245_t50-1450ns_z0-10mm_cmapjet0-15e16_CUT_rotated.png'
    pluto_nframes =list(range(0,256,5))
    l_cap = 1e-2
    extent_measure = [50., 1450., -l_cap/2*1e2, +l_cap/2*1e2]
    t_max = 1250.  # ns
    # This is just for the title label (initial mass density imposed in simulation)
    rho0 = {'mean': 6e-7,
            'max': 2.9e-7}
elif measure_choice=='90A-3cm':
    sim = {'mean': '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-7-I90-3.2cmL-1mmD',
           'max' : '/home/ema/simulazioni/sims_pluto/perTesi/rho3.7e-7-I90-3.2cmL-1mmD-r60-NTOT8-diffRecPeriod8-fast'}
    measurement = '/home/ema/Dottorato/dati_sperimentali_e_calcoli/misure_ne_da_confrontare/misure_capillare_3cmL-1mmD/stessa_scala_temporale/densità_elettronica/ne_I90_t150-1250ns_z0-30mm_cmapjet0-15e16_CUT_rotated.png'
    # ---
    # The frames of pluto which I want to see (it must be a list of integers)
    pluto_nframes =list(range(0,251,5))
    # Capillary radius (m)
    l_cap = 3e-2
    extent_measure = [150., 1250., -l_cap/2*1e2, +l_cap/2*1e2]  # ns, ns, cm, cm
    t_max = 1250.  # ns
    # This is just for the title label (initial mass density imposed in simulation)
    rho0 = {'mean': 8e-7,
            'max': 3.7e-7}
else:
    raise ValueError('Unrecognized measure choice')

average_ne = ('mean', 'max')
#%% Get ne, z-profile
ne_avg_r = {}; t = {}; z={}; z_cc={}; ne_avg_r={}; ne_avg_r_cc={}; tt={}; zz={}
for kind in average_ne:
    z_cc[kind], z[kind], ne_avg_r[kind], t[kind] = map(np.array, apl.ne_avg_over_r(sim[kind], pluto_nframes, kind,
                                                                                   ret_z_cell_borders=True))
    # I average ne over the time cells, so that I get a quantity that is cell centered
    # also on the time grid (before it was vertex centered on the time grid and cell centered
    # on the z grid)
    ne_avg_r_cc[kind] = 0.5*(ne_avg_r[kind][1:,:]+ne_avg_r[kind][:-1,:])

    # Maybe the code would work even without flipping, but I do so, to make if more robust
    z[kind] = np.concatenate((np.flip(-z[kind], axis=0), z[kind][1:]))
    ne_avg_r_cc[kind] = np.concatenate((np.flip(ne_avg_r_cc[kind], axis=1), ne_avg_r_cc[kind]), axis=1)

    tt[kind],zz[kind] = np.meshgrid(t[kind],z[kind])

#%% Read picture
gs = gridspec.GridSpec(6, 2,
                       width_ratios=[30, 2],
                       height_ratios=[1, 1, 1, 1, 1, 1]
                       )
fig = plt.figure(figsize=(5.5,5.5))
ax_sim = [[],[]]
ax_sim[0]  = plt.subplot(gs[2:4,0])
ax_sim[1]  = plt.subplot(gs[4:,0], sharex = ax_sim[0])
ax_meas = plt.subplot(gs[:2,0], sharex = ax_sim[0])
ax_cb = plt.subplot(gs[1:-1,1])

filippi = mpimg.imread(measurement)
# fig, ax = plt.subplots(ncols=2)
ax_meas.imshow(filippi,
               extent=extent_measure,  # extent=[left,right,bottom,top],
               aspect='auto',
               origin='upper'
               )
mp = {}
yticks = [extent_measure[2], 0., extent_measure[3]]
for ii in range(len(average_ne)):
    mp[ii] = ax_sim[ii].pcolormesh(tt[average_ne[ii]], zz[average_ne[ii]], ne_avg_r_cc[average_ne[ii]].T,
                                       cmap='jet',
                                       vmax=1.5e17,
                                       vmin=0.0)

    # ax_sim.invert_yaxis()
    ax_sim[ii].set_ylim([extent_measure[2], extent_measure[3]])
    ax_sim[ii].set_xlim([extent_measure[0], min(extent_measure[1],t_max)])
    ax_sim[ii].set_yticks(yticks)
ax_meas.set_yticks(yticks)

# ax_sim.set_xticks()
# ax_meas.set_xticks(ax_sim.get_xticks())
# ax_meas.set_xticklabels([])

ax_sim[-1].set_xlabel('Time (ns)')
for ii in range(len(ax_sim)):
    ax_sim[ii].set_ylabel('z (cm)')
ax_meas.set_ylabel('z (cm)')

ax_meas.set_title('Measurement')
for ii in range(len(average_ne)):
    title = 'Simulation, ' + r'$\mathrm{' + average_ne[ii] + r'}_r (n_e) $, ' + \
            r'$\rho_0= {:g}'.format(rho0[average_ne[ii]]/1e-7) + r'\cdot 10^{-7}$ g/cm³'
    ax_sim[ii].set_title(title)

cb = fig.colorbar(mp[0], cax = ax_cb, label='Electron density ($\mathrm{cm}^{-3}$)')
cb.set_ticks([0.0, 0.5e17, 1.0e17, 1.5e17])
fig.tight_layout()
