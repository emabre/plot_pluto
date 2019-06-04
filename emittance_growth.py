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
paper_emulate = 'Pompili2017'
# paper_emulate = 'Pompili2018'

# #sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
# ---

pluto_nframes = list(range(0,150,5))  # list(range(0,301,10))
time_unit_pluto = 1e-9  # unit time in pluto's simulation (in s)

# ----- Beam -----
# Normalized emittance (m*rad)
if paper_emulate == 'Pompili2018':
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
    l_cap = 3.2e-2  # m
    r_cap = 0.5e-3  # m
    sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD-r60'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-7-I235-3.2cmL-1mmD'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-8-I235-3.2cmL-1mmD'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/600mbar-I235-3.2cmL-1mmD'
    Dz = 20e-2  # meters

elif paper_emulate == 'Pompili2017':
    emitt_Nx = 1.e-6
    emitt_Ny = emitt_Nx
    energy_MeV = 126
    sigma_x = 130.e-6
    sigma_y = sigma_x
    # d_sigma_x = 0
    d_sigma_x = -(130.-110.)/20.*1.e-4  # Circa... (vedi Fig 6b. sigma(z=9cm)=200um, sigma(z=11cm)=150um -> sigma'=Dsigma/Dz=25*10^-4)
    d_sigma_y = d_sigma_x
    l_cap = 3.2e-2  # m
    r_cap = 0.5e-3  # m
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-8-I90-3.2cmL-1mmD'
    sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I90-3.2cmL-1mmD-r60-NTOT8'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-7-I90-3.2cmL-1mmD'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-6-I90-3.2cmL-1mmD'
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
idx_part_outside_cap = (x**2+y**2 > r_cap**2)
print('{} of {} beam particles are ouside capillary, I remove them.'.format(np.sum(idx_part_outside_cap),
                                                                            Npart))
# x = np.delete(x, idx_part_outside_cap)
# y = np.delete(y, idx_part_outside_cap)

#%% Particles pass in APL, Test case, ideal (no aberration, uniform k)
I = 70.  # Ampere
k_test = cst.mu_0/(2*np.pi) * (cst.e/(cst.m_e*cst.c)) * I/(gamma*r_cap**2)
K_test = k_test*l_cap

# Particles direction change (thin lens approx)
Dxp_test = - K_test*x

xp_new_test = xp+Dxp_test

emitt_x_new_test, sigma_x_new_test, sigma_xp_new_test, cov_xxp_new_test   = apl.emittance(x, xp)
emitt_Nx_test = gamma*emitt_x
emitt_Nx_new_test = gamma*emitt_x_new_test

#%% Particles pass in real APL
times, r_c, g_real, Dg_real = apl.g_Dg_time_evol(sim, pluto_nframes, r_cap, l_cap)

times = times*time_unit_pluto

g_real_interp = np.zeros((len(x), g_real.shape[1]))
for tt in range(g_real.shape[1]):
    g_real_interp[:,tt] = np.interp(np.sqrt(x**2+y**2),
                                    np.concatenate((np.flip(-r_c[1:], axis=0), r_c)),
                                    np.concatenate((np.flip(g_real[1:,tt], axis=0), g_real[:,tt])))

k = cst.e/(cst.m_e*cst.c*gamma) * g_real_interp
K = k*l_cap
Dxp = np.zeros(K.shape)
xp_new = np.zeros(K.shape)
for tt in range(K.shape[1]):
    Dxp[:,tt] = - K[:,tt]*x
    xp_new[:,tt] = xp + Dxp[:,tt]

# New emittance after lens
emitt_x_new = np.zeros(xp_new.shape[1])
for tt in range(K.shape[1]):
    emitt_x_new[tt] = apl.emittance(x, xp_new[:,tt])[0]

# New spot after drift Dz following lens
sigma_x_new = np.zeros(xp_new.shape[1])
x_new = np.zeros(K.shape)
for tt in range(K.shape[1]):
    x_new[:,tt] = x + Dz*(xp_new[:,tt])
    sigma_x_new[tt] = np.std(x_new[:,tt])

# Get current set in simulation
t, I = ut.get_currtab(sim)
#I_apl = np.interp(times, t, I)

emitt_Nx_new = emitt_x_new*gamma

#%% Plot
#plt.close('all')

# Emittance
fig, ax = plt.subplots()
ax.plot(t*1e9, I, '-', color='k', label='Current')
ax.set_ylim(bottom=0.)
ax_emitt = ax.twinx()
ax_emitt.plot(times*1e9, emitt_Nx_new*1e6, 'o-', color='b', label='Emitt.')
ax_emitt.axhline(y=emitt_Nx*1e6, linestyle='--', color='b', label='Emitt. no plasma')
ax_emitt.set_ylabel('Emittance (mm mrad)')
ax_emitt.set_ylim(bottom=0., top=15.)
fig.legend()
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Current (A)')
title = os.path.basename(sim) + "\nσ={:.3g}μm, σ'={:.3g}, ε={:.3g} mm mrad".format(1e6*sigma_x,
                                                                                     d_sigma_x,
                                                                                     emitt_Nx)
title += "Lc={:.3g}cm, Rc={:.3g}mm".format(1e2*l_cap,
                                                 1e3*r_cap)
fig.suptitle(title, color='r')
# ax.set_title(title)
#ax.legend()
plt.tight_layout()

# Spot
fig, ax = plt.subplots()
ax.plot(t*1e9, I, '-', color='k', label='Current')
ax.set_ylim(bottom=0.)
ax_spot = ax.twinx()
ax_spot.plot(times*1e9, sigma_x_new*1e6, 'o-', color='b', label='Spot rms')
# ax_spot.axhline(y=emitt_Nx*1e6, linestyle='--', color='b', label='Emitt. no plasma')
ax_spot.set_ylabel('Spot (mm mrad)')
# ax_spot.set_ylim(bottom=0., top=15.)
fig.legend()
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Current (A)')
title = os.path.basename(sim) + "\nσ={:.3g}μm, σ'={:.3g}, ε={:.3g} mm mrad".format(1e6*sigma_x,
                                                                                     d_sigma_x,
                                                                                     emitt_Nx)
title += "Lc={:.3g}cm, Rc={:.3g}mm".format(1e2*l_cap,
                                           1e3*r_cap)
fig.suptitle(title, color='r')
# ax.set_title(title)
#ax.legend()
plt.tight_layout()

# Trace space
fig, ax = plt.subplots(nrows=2)
ax[0].scatter(x,xp)
ax[1].scatter(x_new,xp_new)
