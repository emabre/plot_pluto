import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import pluto_read_frm as prf

plt.close("all")
importlib.reload(prf)

# <codecell>
# Options
# all_sims = [
#             '/home/ema/simulazioni/sims_pluto/I90/newtransp-rho20',
#             '/home/ema/simulazioni/sims_pluto/I90/asEAAC2017_r22/'
#             ]
all_sims = [
            '/home/ema/simulazioni/sims_pluto/I90/newtransp-rho20',
            ]
legends = [
           'rho20',
           ]
pluto_nframes = [80,160, 197]
# z position of z-const lines (in cm)
z_lines = np.linspace(1e-9,1.8,15)
# Capillary radius
r_cap = 0.5e-3
# Capillary length, half of the real one
l_cap = 1.5e-2

show_legend = True

# <codecell>
# Manipulate the input
# If I set only one simulation, with more frames, than I see all the frames for the same simulation
# "act" stands for "actual"
if len(all_sims)==1 and len(pluto_nframes)>1:
    all_sims_act = all_sims*len(pluto_nframes)
    if len(legends)==1 and show_legend:
        legends_act = legends*len(pluto_nframes)
    elif len(legends)!=len(pluto_nframes) and show_legend:
        raise ValueError('legends should be either in same amount as pluto_nframes or as sims_all, or set show_legend=False')
    else:
        legends_act = legends.copy()

# Load the data
ne_sims = []
ne_avg_sims = []
r_sims = []
z_sims = []
cap = []
times = []
for ii in range(len(all_sims_act)):
    pluto_dir = os.path.join(os.path.expandvars(all_sims_act[ii]),'out')
    q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                # time=125.0,
                                                nframe=pluto_nframes[ii]
                                                )
    times.append(t)
    # Convert r and z to cm
    r /= 1e3
    z /= 1e3
    r_sims.append(r)
    z_sims.append(z)
    ne_sims.append(q["ne"])

    # Build capillary shape (False where there is wall, True elsewere)
    # u_cap = cap_shape(r, z, r_cap, l_cap)
    cap.append(q['interBound']==0e0)

    ne_avg_r = []
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
        print("area_r={}".format(area_r))
        print("integ={}".format(integ))

    ne_avg_sims.append(np.array(ne_avg_r))

# <codecell>
# Plots
# Average ne on fixed z positions
fig_avg, ax_avg = plt.subplots()
for ii in range(len(all_sims_act)):
    ax_avg.plot(z_lines, ne_avg_sims[ii], '.-', label=legends_act[ii]+',t={}'.format(times[ii]))
if show_legend:
    ax_avg.legend()

# Colored map of ne per each simulation
for ii in range(len(all_sims_act)):
    fig_ne, ax_ne = plt.subplots()
    ne_map = ax_ne.pcolormesh(z_sims[ii], r_sims[ii], ne_sims[ii].transpose())
    for line in z_lines:
        ax_ne.axvline(x=line, linestyle='--', color='r')
    capwall = np.ma.masked_array(cap[ii], cap[ii])
    cap_map = ax_ne.pcolormesh(z_sims[ii], r_sims[ii],
                               capwall.transpose(),
                               cmap='Greys_r')
    ax_ne.set_title(legends_act[ii]+',t={}'.format(times[ii]))
    fig_ne.colorbar(ne_map)
