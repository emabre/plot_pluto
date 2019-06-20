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

#plt.close("all")
importlib.reload(prf)
importlib.reload(ut)
importlib.reload(apl)

# electron rest energy in MeV
me_MeV = 0.511

#%% Setting
paper_emulate = 'Pompili2017beam-500A-720A'
# paper_emulate = 'Pompili2018'

# #sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
# ---

pluto_nframes = [5,20,24]  # list(range(0,301,10))
time_unit_pluto = 1e-9  # unit time in pluto's simulation (in s)

# ----- Beam -----
# Normalized emittance (m*rad)
if paper_emulate == 'Pompili2017beam-500A-720A':
    emitt_Nx = 0.8e-6
    emitt_Ny = 0.5e-6
    energy_MeV = 127
    # NB: l'aumento di emitt cambia molto al variare di d_sigma_x
    #sigma_x = 100.e-6
    sigma_x = 110.e-6
    sigma_y = sigma_x
    d_sigma_x = -(113.-105.)/25.*1.e-4
    d_sigma_y = d_sigma_x
    # NB: l'aumento di emitt cambia poco al variare di d_sigma_x (varia anche se decommento qualche riga qui sotto)
    #d_sigma_x -= d_sigma_x*0.5
    #d_sigma_x+-= d_sigma_x*0.5
    #d_sigma_x = 0.0
    l_cap = 1.2e-2  # m
    r_cap = (0.5e-3, 0.6e-3)  # m
    sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD-NEWGRID',
           '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD-NEWGRID']
    Dz = 20e-2  # meters
else :
    raise ValueError('Wrong choice for paper to emulate')

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
r_c = [[],[]]; times = [[],[]];
g_real = [[],[]]; Dg_real = [[],[]]
for ss in range(len(sim)):
    times[ss], r_c[ss], g_real[ss], Dg_real[ss] = apl.g_Dg_time_evol(sim[ss], pluto_nframes, r_cap[ss], l_cap)
    times[ss] = times[ss]*time_unit_pluto

k = []; K = []
g_real_interp = [np.zeros((len(x), g_real[0].shape[1])),
                 np.zeros((len(x), g_real[1].shape[1]))]
xp_new = []; Dxp = []
emitt_x_new = []; emitt_Nx_new = []
sigma_x_new = []; x_new = []
for ss in range(len(sim)):
    for tt in range(g_real[ss].shape[1]):
        g_real_interp[ss][:,tt] = np.interp(np.sqrt(x**2+y**2),
                                            np.concatenate((np.flip(-r_c[ss][1:], axis=0), r_c[ss])),
                                            np.concatenate((np.flip(g_real[ss][1:,tt], axis=0), g_real[ss][:,tt])))

    k.append(cst.e/(cst.m_e*cst.c*gamma) * g_real_interp[ss])
    K.append(k[ss]*l_cap)  # l_cap was used for averaging B (in apl.g_Dg_time_evol) so this is ok!
    Dxp.append(np.zeros(K[ss].shape))
    xp_new.append(np.zeros(K[ss].shape))
    for tt in range(K[ss].shape[1]):
        Dxp[ss][:,tt] = - K[ss][:,tt]*x
        xp_new[ss][:,tt] = xp + Dxp[ss][:,tt]

    # New emittance after lens
    emitt_x_new.append(np.zeros(xp_new[ss].shape[1]))
    for tt in range(K[ss].shape[1]):
        emitt_x_new[ss][tt] = apl.emittance(x, xp_new[ss][:,tt])[0]
    emitt_Nx_new.append(emitt_x_new[ss]*gamma)

    # New spot after drift Dz following lens
    sigma_x_new.append(np.zeros(xp_new[ss].shape[1]))
    x_new.append(np.zeros(K[ss].shape))
    for tt in range(K[ss].shape[1]):
        x_new[ss][:,tt] = x + Dz*(xp_new[ss][:,tt])
        sigma_x_new[ss][tt] = np.std(x_new[ss][:,tt])

# Get current set in simulation
t = [[] for ss in range(len(sim))]
I = [[] for ss in range(len(sim))]
for ss in range(len(sim)):
    t[ss], I[ss] = ut.get_currtab(sim[ss])

#%% Plots
# Trace space
fig, ax = plt.subplots(nrows=2)
ax[0].scatter(x,xp)
ax[1].scatter(x_new,xp_new)

#%% Plot for thesis
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Emittance
fig_I_th, ax_em_th = plt.subplots(figsize=(4.5,3.))

emitt_sim, = ax_em_th.plot(times*1e9, emitt_Nx_new*1e6, color='purple',
                          lw=2,
                          # label='$\epsilon_N$ simulated',
                          zorder=10)

emitt_base = ax_em_th.axhline(y=emitt_Nx*1e6, linestyle='--', lw=2,
                               color='purple',
                               # label='$\epsilon_N$ no plasma',
                               zorder=12)

ax_I_th = ax_em_th.twinx()
curr, = ax_I_th.plot(t*1e9, I, '-', lw=3, c='darkgray', zorder=0)
ax_I_th.set_ylim(bottom=0., top=100.)
ax_I_th.set_zorder(ax_em_th.get_zorder()-1)
ax_em_th.patch.set_visible(False)

ax_em_th.set_ylabel('Emittance (mm mrad)')
ax_em_th.set_ylim(bottom=0., top=14.)
ax_em_th.set_xlim([0.,1200])

# ax_em_th.legend(loc=1)
ax_em_th.legend([curr, emitt_base, emitt_sim],
                ['Current',
                 '$\epsilon_N$ no plasma',
                 '$\epsilon_N$ simulated'])

ax_em_th.set_xlabel('Time (ns)')
ax_I_th.set_ylabel('Current (A)')

plt.tight_layout()

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
# Spot
fig_I_th, ax_sp_th = plt.subplots(figsize=(4.5,3.))

spot_sim, = ax_sp_th.plot(times*1e9, sigma_x_new*1e6,
                          color='purple',
                          lw=2,
                          # label='$\epsilon_N$ simulated',
                          zorder=10)

ax_I_th_sp = ax_sp_th.twinx()
curr, = ax_I_th_sp.plot(t*1e9, I, '-', lw=3, c='darkgray', zorder=0)
ax_I_th_sp.set_ylim(bottom=0., top=100.)
ax_I_th_sp.set_zorder(ax_sp_th.get_zorder()-1)
ax_sp_th.patch.set_visible(False)

ax_sp_th.set_ylabel('Spot rms (μm)')
ax_sp_th.set_xlim([0.,1200])

# ax_sp_th.legend(loc=1)
ax_sp_th.legend([curr, spot_sim],
                ['Current',
                 '$\sigma$ simulated'])

ax_em_th.set_xlabel('Time (ns)')
ax_I_th_sp.set_ylabel('Current (A)')

plt.tight_layout()
