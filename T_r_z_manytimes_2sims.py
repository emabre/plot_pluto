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
pluto_nframes = [5, 10, 20, 80]
# Capillary radius
r_cap = (0.5e-3, 0.6e-3)
# Capillary length, half of the real one
l_cap = 0.6e-2
# Simulated with pluto
sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD-NEWGRID',
        '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD-NEWGRID'
        ]
z_pos = 0.05e-2  # meters
r_pos = 0.00001e-2  # meters
rmax = 0.6e-3 # meters
zmax = 0.8e-2 # meters
name = ('(a)', '(b)')

# <codecell> Load the data
Tr = [[] for ii in range(len(sim))]; Tz = [[] for ii in range(len(sim))]
r_cc = [[] for ii in range(len(sim))]; z_cc = [[] for ii in range(len(sim))]
times = np.zeros((len(sim), len(pluto_nframes)))
for ss in range(len(sim)):
    for ff in range(len(pluto_nframes)):
        pluto_dir = os.path.join(os.path.expandvars(sim[ss]), 'out')
        q, r_temp, z_temp, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                        # time=125.0,
                                                        nframe=pluto_nframes[ff])
        times[ss,ff] = t
        # Convert r and z to m, and compute coordinates of grid centers
        r_cc[ss].append(0.5 * (r_temp[1:]+r_temp[:-1]) / 1e5)
        z_cc[ss].append(0.5 * (z_temp[1:]+z_temp[:-1]) / 1e5)
        Tr[ss].append(q["T"][np.argmin(np.abs(z_cc[ss][ff]-z_pos)), :])
        Tz[ss].append(q["T"][:, np.argmin(np.abs(r_cc[ss][ff]-r_pos))])

# <codecell>
# Plots
fig, ax = plt.subplots(ncols=len(sim), nrows=2, figsize=(5.5, 4.5))

for ss in range(len(sim)):
    for ff in range(len(pluto_nframes)):
        ax[ss,0].plot(r_cc[ss][ff][r_cc[ss][ff]<=r_cap[ss]]*1e6, Tr[ss][ff][r_cc[ss][ff]<=r_cap[ss]],
                      label='{}ns'.format(times[ss,ff]))
        ax[ss,1].plot(z_cc[ss][ff]*1e2, Tz[ss][ff],
                      label='{:.0f}ns'.format(times[ss,ff]))

        ax[ss,0].set_xlim([0,rmax*1e6])
        ax[ss,1].set_xlim([0,zmax*1e2])

        ax[ss,0].set_title(name[ss]+', z={:.0f}cm'.format(z_pos*1e2))
        ax[ss,1].set_title(name[ss]+', r={:.0f}μm'.format(r_pos*1e6))

        # ax[ss,1].set_yticklabels([])

for ss in range(len(sim)):
    for ii in range(ax.shape[1]):
        ax[ss,ii].axvline(x=r_cap[ss]*1e6, ls='-', c='k')

for ii in range(ax.shape[0]):
    ax[ii,0].set_ylabel('Tempearture (K)')
ax[-1,0].set_xlabel('r (μm)')
ax[-1,1].set_xlabel('z (cm)')

for ii in range(ax.shape[0]):
    for jj in range(ax.shape[1]):
        ax[ii,jj].set_ylim([0.,80000])
        ax[ii,jj].grid()

ax[1,1].legend(framealpha=0.4, handlelength=0.5)
fig.tight_layout()



#
