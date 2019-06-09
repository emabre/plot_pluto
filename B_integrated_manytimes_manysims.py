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
# Simulated with pluto
sims = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD',
       '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD'
       ]
pluto_nframes = [10,20,24]  # 10
# Capillary length, half of the real one, including electrodes
l_cap = 0.6e-2  # m
r_cap = (0.5e-3, 0.6e-3)  # m
dz_cap = 0.1e-2  # m
# Static B as in Bobrova model
B_staticBob = np.loadtxt('/home/ema/myprogr/scarica_a_regime/B_adimensional.dat')

# <codecell> Load the data
B = []
r_sims = []
z_sims = []
cap = []
times = []
B_avg_z = []
for ss in range(len(sims)):
    B.append([])
    r_sims.append([])
    z_sims.append([])
    cap.append([])
    times.append([])
    B_avg_z.append([])
    for ii in range(len(pluto_nframes)):
        pluto_dir = os.path.join(os.path.expandvars(sims[ss]), 'out')
        q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                        # time=125.0,
                                                        nframe=pluto_nframes[ii])
        times[ss].append(t)
        # Convert r and z to m
        r /= 1e5
        z /= 1e5
        r_sims[ss].append(r)
        z_sims[ss].append(z)
        B[ss].append(q["bx3"]*1e-4)  # Convert to Tesla

        # Build capillary shape (False where there is wall, True elsewere)
        # u_cap = cap_shape(r, z, r_cap, l_cap)
        cap[ss].append(q['interBound'] == 0e0)

        # Averaging of B over z (I cannot use trapz, since I have cell-averaged B)
        # B_integ_z = np.sum(((B[ss][ii]*cap[ss][ii].astype(int)).T*np.diff(z)).T, axis=0)
        B_integ_z = np.sum(((B[ss][ii]).T*np.diff(z)).T, axis=0)
        B_avg_z[ss].append(B_integ_z/(l_cap-dz_cap))

# <codecell>
# Plots
# Average ne on fixed z positions
fig_avg, ax_avg = plt.subplots(figsize=(3.5, 2.7))
# Fix artifically the last point of B, otherwise bad looking (due to mesh)
# for ss in range(len(sims)):
#     B_avg_z[ss][0][np.argmin(np.abs(500.-0.5*(r[:-1]+r[1:])*1e6))+1] = 1e-3*B_measured[-1,1]
# plot B from PLUTO

# Build cell centered r
r_cc = []
for ss in range(len(sims)):
    r_cc.append([])
    for ii in range(len(pluto_nframes)):
        r_cc[ss].append(0.5*(r_sims[ss][ii][:-1]+r_sims[ss][ii][1:]))
# Add B=0 in r=0 points for better looking plots
for ss in range(len(sims)):
    for ii in range(len(pluto_nframes)):
        B_avg_z[ss][ii] = np.concatenate(([0.], B_avg_z[ss][ii]), axis=0)
        r_cc[ss][ii] = np.concatenate(([0.], r_cc[ss][ii]), axis=0)

linestyles = ['-','--']
colors = ['g', 'darkred']
for ss in range(len(sims)):
    for ii in range(len(pluto_nframes)):
        ax_avg.plot(r_cc[ss][ii][r_cc[ss][ii]<=r_cap[ss]]*1e6,
                    B_avg_z[ss][ii][r_cc[ss][ii]<=r_cap[ss]]*1e3,
                    linestyle = linestyles[ss],
                    color=colors[ss],
                    # label='Simulated',
                    )
for cc in range(len(r_cap)):
    ax_avg.axvline(x=r_cap[cc]*1e6, color='k', linestyle='-', lw=1.)
# ax_avg.set_title('t={}'.format())
# Plot equilibrium model field
# ax_avg.plot(B_staticBob[:,0]*r_cap*1e6/B_staticBob[-1,0],
#             B_staticBob[:,1]*B_measured[-1,1]/B_staticBob[-1,1],
#             '-',
#             color='darkorange',
#             label='Equil. model')

# ax_avg.legend(framealpha = 0.4)
# ax_avg.set_title('t = {}ns, I = {:g}A'.format(times[ii],
#                                             np.interp(times[ii]*1e-9, *ut.get_currtab(sim))))
ax_avg.set_xlim([0., max(r_cap)*1e6])
ax_avg.set_ylim([0., 300.])
ax_avg.set_ylabel('Magnetic field (longitud. avg.) (mT)')
ax_avg.set_xlabel('Radius (μm)')
ax_avg.grid()
fig_avg.tight_layout()
plt.show()
#ax_zt.set_aspect(1)
