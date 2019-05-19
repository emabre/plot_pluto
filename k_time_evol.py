import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import importlib
import os
import pluto_read_frm as prf
from scipy.constants import mu_0

#plt.close("all")
importlib.reload(prf)
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
pluto_nframes = list(range(300))
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

# <codecell> Load the data
B = []
r_sims = []
z_sims = []
cap = []
times = []
B_avg_z = []
for ii in range(len(pluto_nframes)):
    pluto_dir = os.path.join(os.path.expandvars(sim), 'out')
    q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                    # time=125.0,
                                                    nframe=pluto_nframes[ii])
    times.append(t)
    # Convert r and z to m
    r /= 1e5
    z /= 1e5
    r_sims.append(r)
    z_sims.append(z)
    B.append(q["bx3"]/1e4)  # I convert B to Tesla

    # Build capillary shape (False where there is wall, True elsewere)
    # u_cap = cap_shape(r, z, r_cap, l_cap)
    cap.append(q['interBound'] == 0e0)

    # Averaging of B over z (I cannot use trapz, since I have cell-averaged B)
#    B_integ_z = np.sum(((B[ii]*cap[ii].astype(int)).T*np.diff(z)).T,
#                                       axis=0)
    B_integ_z = np.sum(((B[ii]*cap[ii]).T*np.diff(z)).T,
                                       axis=0)
    B_avg_z.append(2*B_integ_z/l_cap)

#%% Build the focusing strenght k(r)=B/r (vertex centered)
times = np.array(times)
# Vertex centered integral of B
B_avg_z_v = []
for ii in range(len(B_avg_z)):
    B_avg_z_v.append(np.concatenate((np.array([0]),
                                     0.5*(B_avg_z[ii][1:] + B_avg_z[ii][:-1]),
                                     np.array([np.nan]))))

# I radially restrict the B field to the capillary
B_avg_z_v_c = np.transpose(np.array([list(B_avg_z_v[ii][r<=r_cap]) for ii in range(len(B_avg_z_v))]))
r_c = r[r<=r_cap]

# Build k (it is not really k, just B/R)
k = np.zeros(B_avg_z_v_c.shape)
for ii in range(B_avg_z_v_c.shape[1]):
    k[1:,ii] = B_avg_z_v_c[1:,ii]/r_c[1:]
k[0,:] = k[1,:]
# Build delta k
Dk = np.zeros(k.shape)
for ii in range(k.shape[1]):
    Dk[:,ii] = k[0,ii] - k[:,ii]

# <codecell> Plots
# Countour of Delta k
#fig, ax = plt.subplots(nrows=2, sharex=True)
    
gs = gridspec.GridSpec(2, 2,
                       width_ratios=[25, 1],
                       height_ratios=[1, 1]
                       )

fig = plt.figure()
ax_Dk = plt.subplot(gs[0,0])
ax_Dk_color = plt.subplot(gs[0,1])
ax_k = plt.subplot(gs[1,0])

tt, rr = np.meshgrid(times, r_c)
lev = np.linspace(-1.e-10,250.,11)
mp = ax_Dk.contourf(tt, rr, Dk, lev, cmap='hot')

#for ii in range(len(r_c)):
#    ax[0].axhline(y=r_c[ii], c='k')

ax_k.set_xlabel('Time (ns)')
ax_Dk.set_ylabel('r ()')
# ax.set_ylim([0.,r_cap*1e+2])
fig.colorbar(mp, cax=ax_Dk_color)

ax_k.plot(times, k[0,:], '-', label=r'$r\rightarrow 0$')
ax_k.plot(times, k[-1,:], '-', label='$r=R$')  # QUESTO NON E' ANCHE PROP. ALLA CORRENTE, PERCHE' USO UN CAMPO INTEGRATO E DIVISO PER L
ax_k.set_ylabel(r'$\frac{<B(r)>}{L}$')
ax_k.legend()

#ax_I = ax_k.twinx()
#ax_I.set_ylim(0., 2*np.pi*r_cap**2/mu_0*ax_k.get_ylim()[1])
#ax_I.set_ylim(2*np.pi*r_cap**2/mu_0*np.array(ax_k.get_ylim()))

# Set lims to have common lims
ax_k.set_xlim(ax_Dk.get_xlim())
ax_k.set_ylim([0,800])  # bottom=0

fig.tight_layout()

plt.show()
