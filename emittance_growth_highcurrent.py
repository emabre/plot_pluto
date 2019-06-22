# -*- coding: utf-8 -*-

import numpy as np
import scipy.constants as cst
import matplotlib.pyplot as plt
from matplotlib import gridspec
import importlib
import os
import pluto_read_frm as prf
from scipy.constants import mu_0
import utilities as ut
import active_plasma_lens as apl

importlib.reload(prf)
importlib.reload(ut)
importlib.reload(apl)

# plt.close('all')

# electron rest energy in MeV
me_MeV = 0.511

#%% Setting
paper_emulate = 'Pompili2017beam-500A-720A'
# paper_emulate = 'Pompili2018'

# #sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
# ---

pluto_nframes = list(range(0,100,1))  # list(range(0,100,5))
frame_track_drift = 40 # Frame of pluto_nframes at which track particles in drift
time_unit_pluto = 1e-9  # unit time in pluto's simulation (in s)

# ----- Beam -----
# Normalized emittance (m*rad)
if paper_emulate == 'Pompili2017beam-500A-720A':
    emitt_Nx = 1.e-6
    emitt_Ny = emitt_Nx
    energy_MeV = 126
    # NB: l'aumento di emitt cambia molto al variare di d_sigma_x
    #sigma_x = 100.e-6
    sigma_x = 130.e-6
    sigma_y = sigma_x
    d_sigma_x = -(130.-110.)/20.*1.e-4  # Circa... (vedi Fig 6b. sigma(z=9cm)=200um, sigma(z=11cm)=150um -> sigma'=Dsigma/Dz=25*10^-4)
    d_sigma_y = d_sigma_x
    # NB: l'aumento di emitt cambia poco al variare di d_sigma_x (varia anche se decommento qualche riga qui sotto)
    #d_sigma_x -= d_sigma_x*0.5
    #d_sigma_x+-= d_sigma_x*0.5
    #d_sigma_x = 0.0
    l_cap = 1.2e-2  # m
    r_cap = (0.5e-3, 0.6e-3)  # m
    sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD-NEWGRID',
           '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD-NEWGRID']
    names = ['(a)', '(b)']
    Dz = 0.1e-2  # meters
    drift_points = np.linspace(0.,20e-2, 100)
else :
    raise ValueError('Wrong choice for paper to emulate')

if frame_track_drift not in pluto_nframes:
    raise ValueError('Choose a frame_track_drift that is in pluto_nframes')

#%% Computations
# Number of particles in beam
Npart = 10000
gamma = energy_MeV/me_MeV

# Build x distribution ---
emitt_x = emitt_Nx/gamma
x, xp = apl.generate_beam_transverse(sigma_x, d_sigma_x, emitt_x, Npart)
# Build y distribution ---
emitt_y = emitt_Ny/gamma
y, yp = apl.generate_beam_transverse(sigma_y, d_sigma_y, emitt_y, Npart)
# Clean particles outside capillary
idx_part_outside_cap = (x**2+y**2 > min(r_cap)**2)
print('{} of {} beam particles are ouside capillary, I remove them.'.format(np.sum(idx_part_outside_cap),
                                                                            Npart))
# x = np.delete(x, idx_part_outside_cap)
# y = np.delete(y, idx_part_outside_cap)

#%% Particles pass in real APL
r_c = len(sim)*[None]; times = len(sim)*[None];
g_real = len(sim)*[None]; Dg_real = len(sim)*[None];
sigma_x_new = len(sim)*[None]; emitt_Nx_new = len(sim)*[None]
emitt_x_new = len(sim)*[None]; x_new = len(sim)*[None]; xp_new = len(sim)*[None]
y_new = len(sim)*[None]; yp_new = len(sim)*[None]
for ss in range(len(sim)):
    times[ss], r_c[ss], g_real[ss], Dg_real[ss] = apl.g_Dg_time_evol(sim[ss], pluto_nframes, r_cap[ss], l_cap)
    times[ss] = times[ss]*time_unit_pluto

    # Map over time
    # (sigma_x_new[ss],
    #  emitt_x_new[ss],
    #  x_new[ss],
    #  xp_new[ss]) = zip(*(map(lambda v: apl.focus_in_thin_apl(v, r_c[ss], x, xp, y, l_cap, gamma, Dz),
    #                          g_real[ss].T)))
    (sigma_x_new[ss],
     emitt_x_new[ss],
     x_new[ss],
     xp_new[ss],
     y_new[ss],
     yp_new[ss]) = zip(*(map(lambda v: apl.focus_in_thick_apl(v, r_c[ss], x, xp, y, yp, l_cap, gamma, Dz, Nz = 10),
                             g_real[ss].T)))

    emitt_Nx_new[ss] = np.array(emitt_x_new[ss])*gamma
    sigma_x_new[ss] = np.array(sigma_x_new[ss])

# Propagate beam through drift after APL
# x_drift, sigma_x_drift = map(apl.drift_beam(x_new, xp_new, drift_points)
x_drift = len(sim)*[None]; sigma_x_drift= len(sim)*[None]
tt = np.argwhere(frame_track_drift==np.array(pluto_nframes))[0,0]
for ss in range(len(sim)):
    x_drift[ss], sigma_x_drift[ss] = apl.drift_beam(x_new[ss][tt], xp_new[ss][tt], drift_points)


# Get current set in simulation
t = [[] for ss in range(len(sim))]
I = [[] for ss in range(len(sim))]
for ss in range(len(sim)):
    t[ss], I[ss] = ut.get_currtab(sim[ss])

#%% Plot for thesis
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
linestyles = ('--', '-')
# Emittance
fig_I_th, ax_em_th = plt.subplots(figsize=(4.5,3.))
ax_I_th = ax_em_th.twinx()
curr = []; emitt_sim = []
for ss in range(len(sim)):
    emitt_sim.append(ax_em_th.plot(times[ss]*1e9, emitt_Nx_new[ss]*1e6, color='purple',
                              lw=2,
                              ls = linestyles[ss],
                              # label = names[ss] + ', $\epsilon_N$',
                              zorder = 10)[0])
    curr.append(ax_I_th.plot(t[ss]*1e9, I[ss],
                         linestyles[ss],
                         lw=3,
                         c='darkgray',
                         # label = names[ss]+', current',
                         zorder=0)[0])

ax_I_th.set_ylim(bottom=0., top=800.)
ax_I_th.set_zorder(ax_em_th.get_zorder()-1)
ax_em_th.patch.set_visible(False)

ax_em_th.set_ylabel('Emittance (mm mrad)')
ax_em_th.set_ylim(bottom=0., top=5.)
ax_em_th.set_xlim([0.,500])

# ax_em_th.legend(curr + emitt_sim + [emitt_base],
#                 [names[ss]+',current' for ss in range(len(sim))] + \
#                 [names[ss]+',emittance' for ss in range(len(sim))] + \
#                 ['no plasma, emittance'])
ax_em_th.legend(curr + emitt_sim,
                [names[ss]+',current' for ss in range(len(sim))] + \
                [names[ss]+',emittance' for ss in range(len(sim))],
                loc = 4)

ax_em_th.set_xlabel('Time (ns)')
ax_I_th.set_ylabel('Current (A)')

plt.tight_layout()

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Spot as function of z, fixed timing
fig_spz_th, ax_spz_th = plt.subplots(figsize=(4.5,3.))
# Plot spot
for ss in range(len(sim)):
    spot_sim, = ax_spz_th.plot(drift_points*1e2, sigma_x_drift[ss]*1e6,
                               # color='purple',
                               lw=2,
                               # label='$\epsilon_N$ simulated',
                               zorder=10,
                               label=names[ss])

ax_spz_th.set_ylabel('Spot rms (μm)')
ax_spz_th.set_xlim([0.,30.])
ax_spz_th.legend()
ax_em_th.set_xlabel('z (cm)')

fig_spz_th.tight_layout()


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Spot
fig_I_th, ax_sp_th = plt.subplots(figsize=(4.5,3.))
ax_I_th_sp = ax_sp_th.twinx()
# Plot spot
for ss in range(len(sim)):
    spot_sim, = ax_sp_th.plot(times[ss]*1e9, sigma_x_new[ss]*1e6,
                              color='purple',
                              lw=2,
                              # label='$\epsilon_N$ simulated',
                              zorder=10)

    curr, = ax_I_th_sp.plot(t[ss]*1e9, I[ss], '-', lw=3, c='darkgray', zorder=0)
# ax_I_th_sp.set_ylim(bottom=0., top=100.)
ax_I_th_sp.set_zorder(ax_sp_th.get_zorder()-1)
ax_sp_th.patch.set_visible(False)

ax_sp_th.set_ylabel('Spot rms (μm)')
ax_sp_th.set_xlim([0.,500])

# ax_sp_th.legend(loc=1)
ax_sp_th.legend([curr, spot_sim],
                ['Current',
                 '$\sigma$ simulated'])

ax_em_th.set_xlabel('Time (ns)')
ax_I_th_sp.set_ylabel('Current (A)')

plt.tight_layout()
