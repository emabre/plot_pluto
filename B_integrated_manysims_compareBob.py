import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import pluto_read_frm as prf
import utilities as ut
from scipy import constants as cst

#plt.close("all")
importlib.reload(prf)
# plt.ion()

# <codecell>
# Settings
# Simulated with pluto
sims = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD',
       '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD'
       ]
pluto_nframes = [30]  # 10
# Capillary length, half of the real one, including electrodes
l_cap = 0.6e-2  # m
r_cap = (0.5e-3, 0.6e-3)  # m
dz_cap = 0.1e-2  # m
# Static B as in Bobrova model
B_staticBob = np.loadtxt('/home/ema/myprogr/scarica_a_regime/B_adimensional.dat')
plot_staticBob = True

names = ('(a)', '(b)')

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

# Plot equilibrium model field
scale_factors = [1.09, 1.09]
if plot_staticBob:
    color_Bob = 'darkorange'
    linestyle_Bob = '-'
    for ss in range(len(sims)):
        for ii in range(len(pluto_nframes)):
            # Get I by interpolation on the current table
            t_I, I = ut.get_currtab(sims[ss])
            t0_ns = t_I[np.max(np.argwhere(t_I*1e9-times[ss]<0))]*1e9
            I0_A = I[np.max(np.argwhere(t_I*1e9-times[ss]<0))]
            t1_ns = t_I[np.min(np.argwhere(t_I*1e9-times[ss]>0))]*1e9
            I1_A = I[np.min(np.argwhere(t_I*1e9-times[ss]>0))]
            Inow_A = I0_A + (I1_A-I0_A)/(t1_ns-t0_ns)*(times[ss][ii]-t0_ns)
            Bwall_T = cst.mu_0*Inow_A/(2*np.pi*r_cap[ss])

            ax_avg.plot(B_staticBob[:,0]*r_cap[ss]*1e6/B_staticBob[-1,0],
                        B_staticBob[:,1]*Bwall_T*1e3/B_staticBob[-1,1] *scale_factors[ss],
                        linestyle = linestyle_Bob,
                        linewidth = 2.3,
                        # marker = '^',
                        # markersize = 2.,
                        # markevery = 5,
                        color=color_Bob,
                        # label='Equil. model',
                        )


linestyles = ['-','--']
# colors = [['red', 'darkred'], ['royalblue','navy']]

for ss in range(len(sims)):
    for ii in range(len(pluto_nframes)):
        ax_avg.plot(r_cc[ss][ii][r_cc[ss][ii]<=r_cap[ss]]*1e6,
                    B_avg_z[ss][ii][r_cc[ss][ii]<=r_cap[ss]]*1e3,
                    linestyle = linestyles[ss],
                    linewidth = 0.9,
                    color='k',
                    # label='t={:.0f}ns, R={:.0f}μm'.format(times[ss][ii], r_cap[ss]*1e6),
                    )
# Set legend
import matplotlib.lines as mlines
handles=[]
for ss in range(len(sims)):
    handles.append(mlines.Line2D([], [], color='k', linestyle=linestyles[ss],  # , marker='*', markersize=15,
                   label=names[ss]))
if plot_staticBob:
    handles.append(mlines.Line2D([], [], color=color_Bob, linestyle=linestyle_Bob,  # , marker='*', markersize=15,
                   label='Equil. model'))
ax_avg.legend(handles=handles)

# Set title
ax_avg.set_title("t={:.0f}ns".format(times[0][ii]))
# Set line to separate end of one capillary
for cc in range(len(r_cap)):
    ax_avg.axvline(x=r_cap[cc]*1e6, color='k', linestyle='-', lw=1.)

# fig_avg.text(0.65, 0.8, "t={:.0f}ns,".format(times[ii][0]), ha="center", va="bottom", size="medium",color="navy")
ax_avg.set_xlim([0., max(r_cap)*1e6])
ax_avg.set_ylim([0., 300.])
ax_avg.set_ylabel('Magnetic field (longitud. avg.) (mT)')
ax_avg.set_xlabel('r (μm)')
ax_avg.grid()
fig_avg.tight_layout()
plt.show()
#ax_zt.set_aspect(1)
