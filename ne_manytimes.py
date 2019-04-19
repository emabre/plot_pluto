import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import pluto_read_frm as prf

# plt.close("all")
importlib.reload(prf)

# <codecell>
# Options
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa'
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-rhounif-I90-3cm'
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-rhounif'
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa'
# sim = '/home/ema/simulazioni/sims_pluto/dens_real/1e4Pa-06012019'

# legend = '1.e5Pa'
legend = os.path.split(sim)[1]
# legend = '1.3e5Pa-rhounif'

# How to transversely average ne: 'max' (max value), 'integral' (integral average)
average_ne = 'integral'

plot_ne_map_each_frame = False

# The frames of pluto which I want to see (it must be a list of integers, with
# dimension : len(all_sims)*(number of frames you want to watch in every simulation))
# fastest varying index: frames for the same simulation, slower running index: simulation
#pluto_nframes = [80, 160, 200]
pluto_nframes = [20*ii for ii in range(1,15)]
# pluto_nframes = [80]
# z position of z-const lines (in cm)
# Z lines settings, z lines always start from 0
N_z_lines = 30
z_lines_end = 1.5

show_legend = True

reflect_lowz = True
zlim_plot = 0.75

# <codecell>
# Manipulate the input

z_lines, dz = np.linspace(1e-12, z_lines_end, N_z_lines, retstep=True)

# Load the data
ne_sims = []
ne_avg_sims = []
r_sims = []
z_sims = []
cap = []
times = []
for ii in range(len(pluto_nframes)):
    pluto_dir = os.path.join(os.path.expandvars(sim),'out')
    q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                    # time=125.0,
                                                    nframe=pluto_nframes[ii])
    times.append(t)
    # Convert r and z to cm
    r /= 1e3
    z /= 1e3
    r_sims.append(r)
    z_sims.append(z)
    ne_sims.append(q["ne"])

    ne_avg_r = []
    if average_ne == 'integral':
        # Build capillary shape (False where there is wall, True elsewere)
        cap.append(q['interBound']==0e0)
        areas = []
        for jj in range(len(z_lines)):
            # Integration on line z=const
            # z_lines will be between z[idx_z] and z[idx_z+1]
            idx_z = np.argmax(z[z_lines[jj]>=z])
            # Just to debug
            # print("z[idx_z]={};\tz_line={};\tz[idx_z+1]={}".format(z[idx_z],
            #                                                      z_lines[jj],
            #                                                      z[idx_z+1]))
            integ = np.sum(np.pi * q["ne"][idx_z,:] * (r[1:]**2 - r[:-1]**2) * cap[ii][idx_z,:])
            area_r = np.sum(np.pi * (r[1:]**2 - r[:-1]**2) * cap[ii][idx_z,:])
            areas.append(area_r)
            ne_avg_r.append(integ/area_r)
            # print("area_r={}".format(area_r))
            # print("integ={}".format(integ))
    elif average_ne == 'max':
        for jj in range(len(z_lines)):
            idx_z = np.argmax(z[z_lines[jj]>=z])
            ne_avg_r.append(q["ne"][idx_z,:].max())
    else:
        raise ValueError('Wrong choice for average_ne')

    ne_avg_sims.append(np.array(ne_avg_r))

if reflect_lowz:
    # Mayve the code would work even without flipping, but I do so, to make if more robust
    z_lines = np.append(np.flip(-z_lines, axis=0), z_lines)
    for ii in range(len(pluto_nframes)):
        ne_avg_sims[ii] = np.append(np.flip(ne_avg_sims[ii], axis=0), ne_avg_sims[ii])

#z_lines = z_lines+0.5
#z_lines*=10

# <codecell>
# Plots
# Average ne on fixed z positions
fig_avg, ax_avg = plt.subplots()

for ii in range(len(pluto_nframes)):
    ax_avg.plot(z_lines, ne_avg_sims[ii], '.-', label=legend+',t={}'.format(times[ii]))
    ax_avg.grid()
if show_legend:
    ax_avg.legend()
ax_avg.set_ylabel('$n_e / \mathrm{cm}^{-3}$ (transverse average)')
ax_avg.set_xlabel('$z / \mathrm{cm}$')
ax_avg.set_xlim([-zlim_plot, zlim_plot])
#ax_avg.set_xlim([0.0, 10.0])


if plot_ne_map_each_frame:
    # Colored map of ne per each simulation frame
    for ii in range(len(pluto_nframes)):
        fig_ne, ax_ne = plt.subplots()
        ne_map = ax_ne.pcolormesh(z_sims[ii], r_sims[ii], ne_sims[ii].transpose())
        for line in z_lines:
            ax_ne.axvline(x=line, linestyle='--', color='r')
        capwall = np.ma.masked_array(cap[ii], cap[ii])
        cap_map = ax_ne.pcolormesh(z_sims[ii], r_sims[ii],
                                   capwall.transpose(),
                                   cmap='Greys_r')
        ax_ne.set_title(legend+',t={}'.format(times[ii]))
        fig_ne.colorbar(ne_map)

# Colored map of ne, where x is longitudinal coordinate and y is time
# (like in Francesco's thesis)
ne_zt = np.stack(ne_avg_sims, axis=1).transpose()
fig_zt, ax_zt = plt.subplots()
ne_map = ax_zt.imshow(ne_zt, origin='upper',
                     extent=[z_lines[0]-dz*0.5, z_lines[-1]+dz*0.5,
                             times[-1]+0.5*(times[-1]-times[-2]), times[0]-0.5*(times[1]-times[0])],  # extent=[left,right,bottom,top]
                     aspect='auto',
                     cmap='jet',
                     vmax = 1.5e17
                     )
fig_zt.colorbar(ne_map)
ax_zt.set_title(legend + '\ntransv. avg: ' + average_ne)
plt.show()
#ax_zt.set_aspect(1)
