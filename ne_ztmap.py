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
measure_choice = '245A-1cm'  # '90A-3cm' or '245A-1cm'
average_ne = 'max'  # 'max','integral'
# This is just for the title label
rho0 = 3.1e-7   # Initial mass density imposed in simulation

if measure_choice=='245A-1cm':
    # ----
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1e5Pa'
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1e4Pa-06012019/'
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa'
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho1.5e-7-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho3.5e-7-I245-1.2cmL-1mmD-NEWGRID/'

    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho6e-7-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-7-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho1e-6-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho1.5e-6-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho1.75e-6-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2e-6-I245-1.2cmL-1mmD-NEWGRID'

    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho3.3e-7-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho3.1e-7-I245-1.2cmL-1mmD-NEWGRID'
    sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.9e-7-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho4e-6-I245-1.2cmL-1mmD-NEWGRID'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I245-1.2cmL-1mmD-NEWGRID/'
    measurement = '/home/ema/Dottorato/dati_sperimentali_e_calcoli/misure_ne_da_confrontare/misure_capillare_1cmL-1mmD/stessa_scala_temporale/ne_I245_t50-1450ns_z0-10mm_cmapjet0-15e16_CUT_rotated.png'
    pluto_nframes =list(range(0,150,5))
    l_cap = 1e-2
    extent_measure = [50., 1450., -l_cap/2*1e2, +l_cap/2*1e2]
elif measure_choice=='90A-3cm':
    # ---
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-8-I90-3.2cmL-1mmD'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I90-3.2cmL-1mmD-r60-NTOT8'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho3.5e-7-I90-3.2cmL-1mmD-r60-NTOT8-diffRecPeriod8-fast'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/100mbarOK-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/300mbarOK-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
    sim = '/home/ema/simulazioni/sims_pluto/perTesi/200mbarOKselfmade-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho4.5e-7-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-7-I90-3.2cmL-1mmD'
    measurement = '/home/ema/Dottorato/dati_sperimentali_e_calcoli/misure_ne_da_confrontare/misure_capillare_3cmL-1mmD/stessa_scala_temporale/densità_elettronica/ne_I90_t150-1250ns_z0-30mm_cmapjet0-15e16_CUT_rotated.png'
    # ---
    # The frames of pluto which I want to see (it must be a list of integers)
    pluto_nframes =list(range(0,81,5))
    # Capillary radius (m)
    l_cap = 3e-2
    extent_measure = [150., 1250., -l_cap/2*1e2, +l_cap/2*1e2]  # ns, ns, cm, cm
else:
    raise ValueError('Unrecognized measure choice')

#%% Get ne, z-profile
z_cc, z, ne_avg_r, t = apl.ne_avg_over_r(sim, pluto_nframes, average_ne,
                                              ret_z_cell_borders=True)
t = np.array(t)
ne_avg_r = np.array(ne_avg_r)
# I average ne over the time cells, so that I get a quantity that is cell centered
# also on the time grid (before it was vertex centered on the time grid and cell centered
# on the z grid)
ne_avg_r_cc = 0.5*(ne_avg_r[1:,:]+ne_avg_r[:-1,:])

# Maybe the code would work even without flipping, but I do so, to make if more robust
z = np.concatenate((np.flip(-z, axis=0), z[1:]))
ne_avg_r_cc = np.concatenate((np.flip(ne_avg_r_cc, axis=1), ne_avg_r_cc), axis=1)

tt,zz = np.meshgrid(t,z)

#%% Read picture
gs = gridspec.GridSpec(2, 2,
                       width_ratios=[30, 2],
                       height_ratios=[1, 1]
                       )
fig = plt.figure(figsize=(5.5,4.5))
ax_sim  = plt.subplot(gs[1,0])
ax_meas = plt.subplot(gs[0,0])
ax_cb = plt.subplot(gs[:,1])

filippi = mpimg.imread(measurement)
# fig, ax = plt.subplots(ncols=2)
ax_meas.imshow(filippi,
             extent=extent_measure,  # extent=[left,right,bottom,top],
             aspect='auto',
             origin='upper'
             )
mp = ax_sim.pcolormesh(tt, zz, ne_avg_r_cc.T,
                   cmap='jet',
                   vmax=1.5e17,
                   vmin=0.0)

# ax_sim.invert_yaxis()
ax_sim.set_ylim([extent_measure[2], extent_measure[3]])
ax_sim.set_xlim([extent_measure[0], extent_measure[1]])
yticks = [extent_measure[2], 0., extent_measure[3]]
ax_meas.set_yticks(yticks)
ax_sim.set_yticks(yticks)
# ax_sim.set_xticks()
# ax_meas.set_xticks(ax_sim.get_xticks())
# ax_meas.set_xticklabels([])

ax_sim.set_xlabel('Time (s)')
ax_sim.set_ylabel('z (cm)')
ax_meas.set_ylabel('z (cm)')
if average_ne == 'max':
    title = 'Simulation, '+r'$\mathrm{max}_r (n_e) $, '+r'$\rho_0= {:g}'.format(rho0/1e-7)+r'\cdot 10^{-7}$ g/cm³'
    # title = 'Simulation, '+r'$\mathrm{max}_r (n_e) $, '+'300mbar'
elif average_ne == 'integral':
    # title = 'Simulation, '+r'$\langle n_e \rangle_r$, '+'300mbar'
    title = 'Simulation, '+r'$\langle n_e \rangle_r$, '+r'$\rho_0= {:g}'.format(rho0/1e-7)+r'\cdot 10^{-7}$ g/cm³'
ax_sim.set_title(title)
ax_meas.set_title('Measurement')

cb = fig.colorbar(mp, cax = ax_cb, label='Electron density ($\mathrm{cm}^{-3}$)')
cb.set_ticks([0.0, 0.5e17, 1.0e17, 1.5e17])
fig.tight_layout()
