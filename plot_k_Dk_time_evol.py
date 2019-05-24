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
sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'

# legend = '1.e5Pa'
legend = '1.3e5Pa'

plot_ne_map_each_frame = False

# The frames of pluto which I want to see (it must be a list of integers, with
# dimension : len(all_sims)*(number of frames you want to watch in every simulation))
# fastest varying index: frames for the same simulation, slower running index: simulation
# pluto_nframes = [80, 160, 200]
# pluto_nframes = [-10+20*ii for ii in range(1,15)]
# pluto_nframes = [24, 50, 74, 100, 120, 150]
# pluto_nframes = [24, 50, 74, 100, 120, 132, 150]
pluto_nframes = list(range(321))
# z position of z-const lines (in cm)
# Z lines settings, z lines always start from 0
N_z_lines = 30
z_lines_end = 0.5
# Capillary radius
r_cap = 0.5e-3
# Capillary length, half of the real one
l_cap = 0.5e-2

show_legend = False

reflect_lowz = True
zlim_plot = 0.5
ne_lim_plot = 1.5e17

#%% Compute/get k, Dk and I
times, r_c, g, Dg = apl.g_Dg_time_evol(sim, pluto_nframes, r_cap, l_cap)
t_I, I = ut.get_currtab(sim)

# <codecell> Plots
# Countour of Delta k
#fig, ax = plt.subplots(nrows=2, sharex=True)

gs = gridspec.GridSpec(2, 2,
                       width_ratios=[30, 1],
                       height_ratios=[1, 1]
                       )

fig = plt.figure()
ax_Dk = plt.subplot(gs[0,0])
ax_Dk_color = plt.subplot(gs[0,1])
ax_k = plt.subplot(gs[1,0])

tt, rr = np.meshgrid(times, r_c)
lev = np.linspace(-1.e-10,250.,11)
mp = ax_Dk.contourf(tt, rr*1e6, Dk, lev, cmap='hot')

#for ii in range(len(r_c)):
#    ax[0].axhline(y=r_c[ii], c='k')

ax_k.set_xlabel('Time (ns)')
ax_Dk.set_ylabel('$r$ (μm)')
# ax.set_ylim([0.,r_cap*1e+2])
fig.colorbar(mp, cax=ax_Dk_color, label=r'$\langle g(0) \rangle - \langle g(r) \rangle$ (T/m)')

ax_k.plot(times, k[0,:], '-', label=r'$\langle g \rangle (0)$')
ax_k.plot(times, k[-1,:], '-', label=r'$\langle g \rangle (R)$')  # QUESTO NON E' ANCHE PROP. ALLA CORRENTE, PERCHE' USO UN CAMPO INTEGRATO E DIVISO PER L
ax_k.set_ylabel(r'$\langle g \rangle$'+' (T/m)')

ax_I = ax_k.twinx()
ax_I.plot(t_I*1e9, I, '--', color='k', label='current')
ax_I.set_ylabel('Current (A)')
ax_I.set_ylim(0., 250.)
# ax_I.set_ylim(0., 2*np.pi*r_cap**2/mu_0*ax_k.get_ylim()[1])
# ax_I.set_ylim(2*np.pi*r_cap**2/mu_0*np.array(ax_k.get_ylim()))
fig.legend(loc=(0.6,0.31))



# Set lims to have common lims
ax_k.grid()
ax_k.set_xlim(ax_Dk.get_xlim())
ax_k.set_ylim([0,800])  # bottom=0

fig.tight_layout()

plt.show()
