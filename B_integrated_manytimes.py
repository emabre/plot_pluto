import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import pluto_read_frm as prf

#plt.close("all")
importlib.reload(prf)
plt.ion()

# <codecell>
# Options
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1e5Pa'
sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa'

# legend = '1.e5Pa'
legend = '1.3e5Pa'

plot_ne_map_each_frame = False

# The frames of pluto which I want to see (it must be a list of integers, with
# dimension : len(all_sims)*(number of frames you want to watch in every simulation))
# fastest varying index: frames for the same simulation, slower running index: simulation
# pluto_nframes = [80, 160, 200]
# pluto_nframes = [-10+20*ii for ii in range(1,15)]
pluto_nframes = [50, 74, 100, 120, 150, 190]
# z position of z-const lines (in cm)
# Z lines settings, z lines always start from 0
N_z_lines = 30
z_lines_end = 0.5
# Capillary radius
r_cap = 0.5e-3
# Capillary length, half of the real one
l_cap = 0.5e-2

show_legend = True

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
    # Convert r and z to cm
    r /= 1e3
    z /= 1e3
    r_sims.append(r)
    z_sims.append(z)
    B.append(q["bx3"])

    # Build capillary shape (False where there is wall, True elsewere)
    # u_cap = cap_shape(r, z, r_cap, l_cap)
    cap.append(q['interBound'] == 0e0)

    # Averaging of B over z (I cannot use trapz, since I have cell-averaged B)
    B_integ_z = np.sum(((B[ii]*cap[ii].astype(int)).T*np.diff(z)).T, axis=0)
    B_avg_z.append(B_integ_z/(l_cap/1e-2))

# <codecell>
# Plots
# Average ne on fixed z positions
fig_avg, ax_avg = plt.subplots()

for ii in range(len(pluto_nframes)):
    ax_avg.plot(0.5*(r[:-1]+r[1:]), B_avg_z[ii],
                '.-', label=legend+',t={}'.format(times[ii]))
if show_legend:
    ax_avg.legend()
ax_avg.set_xlim([0.,r_cap/1e-2])
ax_avg.set_ylabel('$B / \mathrm{g}$ (longitudinal average)')
ax_avg.set_xlabel('$r / \mathrm{cm}$')
ax_avg.grid()

plt.show()
#ax_zt.set_aspect(1)
