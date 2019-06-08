import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import pluto_read_frm as prf
import utilities as ut

#plt.close("all")
importlib.reload(prf)
# plt.ion()

# <codecell>
# Settings
pluto_nframes = [60, 120]
# Capillary radius
r_cap = 0.5e-3
# Capillary length, half of the real one
l_cap = 1.5e-2
# Simulated with pluto
sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD-Twall2300',
        '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD',
        '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I235-3.2cmL-1mmD-Twall4600'
        ]
z_pos = 0.8e-2  # meters

# <codecell> Load the data
T = [[] for ii in range(len(sim))]
times = np.zeros((len(sim), len(pluto_nframes)))
for ss in range(len(sim)):
    for ff in range(len(pluto_nframes)):
        pluto_dir = os.path.join(os.path.expandvars(sim[ss]), 'out')
        q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                        # time=125.0,
                                                        nframe=pluto_nframes[ff])
        times[ss,ff] = t
        # Convert r and z to m
        r /= 1e5
        z /= 1e5
        T[ss].append(q["T"][np.argmin(np.abs(z-z_pos)), :])

# <codecell>
# Plots
fig_avg, ax = plt.subplots(ncols=len(pluto_nframes), figsize=(4.5, 2.7))

linestyles = ['-','-','-']
markers = ['o','^','s']
marker_colors = ['yellow', 'orange', 'red']
labels = ['$T_w$=2300K', '$T_w$=3400K', '$T_w$=4600K']
markeveries = [slice(2*ss, -1, 2*len(sim)) for ss in range(len(sim))]
for ss in range(len(sim)):
    for ff in range(len(pluto_nframes)):
        ax[ff].plot(0.5*(r[:-1]+r[1:])*1e6, T[ss][ff],
                    linestyle = linestyles[ss],
                    color='k',
                    label = labels[ss],
                    marker = markers[ss],
                    markevery = markeveries[ss],
                    mfc = marker_colors[ss],
                    mec = 'k', # marker_colors[ss],
                    ms = 4.,
                    lw = 0.5
                    )
        ax[ff].set_title('t = {:g}ns, I = {:g}A'.format(times[0][ff],
                                                      np.interp(times[0][ff]*1e-9, *ut.get_currtab(sim[ss]))))
ax[0].legend(framealpha=0.4)
for ss in range(len(sim)):
    for ff in range(len(pluto_nframes)):
        ax[ff].plot(0.5*(r[:-1]+r[1:])*1e6, T[ss][ff],
                    linestyle = '',
                    label = labels[ss],
                    marker = markers[ss],
                    markevery = markeveries[ss],
                    mfc = marker_colors[ss],
                    mec = 'k', # marker_colors[ss],
                    ms = 4.,
                    lw = 0.5
                    )
        ax[ff].set_xlim([0.,r_cap*1e6])
        ax[ff].set_xticks([0., 250., 500.])
        ax[ff].set_yticks([0., 20000, 40000, 60000])
        ax[ff].set_ylim([0.,60000])
        ax[ff].grid()
        ax[ff].set_xlabel('Radius (μm)')
ax[0].set_ylabel('Temperature (K)')
ax[1].set_yticklabels([])
ax[1].set_yticklabels([])
# ax_avg.set_ylim([0., 25.])
# ax_avg.set_ylabel('Magnetic field (longitud. avg.) (mT)')
# ax_avg.set_xlabel('Radius (μm)')
# ax_avg.grid()
fig_avg.tight_layout()
# plt.show()
# ax_zt.set_aspect(1)
