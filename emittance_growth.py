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

#%% Settings

# ----- Beam -----
# Normalized emittance (m*rad)
emitt_N = 1e-6
energy_MeV = 127
me_MeV = 0.511
# Sigma
sigma_x = 110.e-6
# Derivative of sigma w.r.t. z
d_sigma_x = (113.-105.)/25.*1.e-4
# Number of particles
Npart = 10000

# ----- Simulation -----
sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
pluto_nframes = [25, 50, 75, 100, 125, 150, 175, 200, 225, 250, 275, 300]  # list(range(321))
l_cap = 3e-2  # m
r_cap = 0.5e-3  # m
time_unit_pluto = 1e-9  # unit time in pluto's simulation (in s)

#%% Computations
gamma = energy_MeV/me_MeV
cov_xxp = sigma_x * d_sigma_x
emitt = emitt_N/gamma
sigma_xp = np.sqrt((emitt**2 + cov_xxp**2)/sigma_x**2)

mean = [0,0]
cov = [[sigma_x**2, cov_xxp], [cov_xxp, sigma_xp**2]]

# Generate distribution
x, xp = np.random.multivariate_normal(mean, cov, Npart).T

Npart_outside_cap = np.sum(np.abs(x) > r_cap)
print('{} of {} beam particles are ouside capillary, I remove them.'.format(Npart_outside_cap, Npart))

#%% Check if distro is ok
#cov_xxp_test = np.cov(np.stack([x,xp]))[0,1]
#sigma_x_test = x.std()
#sigma_xp_test = xp.std()
#emitt_N_test = gamma*np.sqrt(sigma_x_test**2 * sigma_xp_test**2 -  cov_xxp_test**2)

emitt_test, sigma_x_test, sigma_xp_test, cov_xxp_test   = apl.emittance(x, xp)
emitt_N_test = gamma*emitt_test

print('generated distro with {} partilces'.format(Npart))
print('Emittance (normalized):')
print('{:.5g} mm mrad (required);{:.5g} mm mrad (obtained)'.format(emitt_N*1e6,
                                                             emitt_N_test*1e6))
print('Spot:')
print('{:.5g} μm (required);{:.5g} μm mrad (obtained)'.format(sigma_x*1e6,
                                                             sigma_x_test*1e6))
print("Covariance x,x':")
print('{:.5g} m (required);{:.5g} m mrad (obtained)'.format(cov_xxp,
                                                             cov_xxp_test))

#%% Particles pass in APL, Test case, ideal (no aberration, uniform k)
I = 70.  # Ampere
k_test = cst.mu_0/(2*np.pi) * (cst.e/(cst.m_e*cst.c)) * I/(gamma*r_cap**2)
K_test = k_test*l_cap

# Particles direction change (thin lens approx)
Dxp_test = - K_test*x

xp_new_test = xp+Dxp_test

emitt_new_test, sigma_x_new_test, sigma_xp_new_test, cov_xxp_new_test   = apl.emittance(x, xp)
emitt_N_test = gamma*emitt_test
emitt_N_new_test = gamma*emitt_new_test

#%% Particles pass in real APL
times, r_c, g_real, Dg_real = apl.g_Dg_time_evol(sim, pluto_nframes, r_cap, l_cap)

times = times*time_unit_pluto

g_real_interp = np.zeros((len(x), g_real.shape[1]))
for tt in range(g_real.shape[1]):
    g_real_interp[:,tt] = np.interp(x,
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
emitt_new = np.zeros(xp_new.shape[1])
sigma_x_new = np.zeros(xp_new.shape[1])
sigma_xp_new = np.zeros(xp_new.shape[1])
cov_xxp_new = np.zeros(xp_new.shape[1])
for tt in range(K.shape[1]):
    (emitt_new[tt],
     sigma_x_new[tt],
     sigma_xp_new[tt],
     cov_xxp_new[tt]) = apl.emittance(x, xp_new[:,tt])

t, I = ut.get_currtab(sim)
#I_apl = np.interp(times, t, I)


emitt_N_new = emitt_new*gamma

#%% Plot
plt.close('all')

fig, ax = plt.subplots()

ax.plot(t*1e9, I, '-', color='k', label='Current')
ax_emitt = ax.twinx()
ax_emitt.plot(times*1e9, emitt_N_new*1e6, '.-', color='b', label='Emitt.')
ax_emitt.axhline(y=emitt_N*1e6, linestyle='--', color='b', label='Emitt. no plasma')
ax_emitt.set_ylabel('Emittance (mm mrad)')
fig.legend()
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Current (A)')
#ax.legend()
plt.tight_layout()