import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import pluto_read_frm as prf

plt.close("all")
importlib.reload(prf)

# <codecell>
# Options
all_sims = [
            '/home/ema/simulazioni/sims_pluto/I90/newtransp-rho20',
            '/home/ema/simulazioni/sims_pluto/I90/asEAAC2017_r22/'
            ]
legends = [
           'rho20',
           'rho as EAAC2017'
           ]
pluto_nframes = [110,110]
# z position of z-const lines (in cm)
z_lines = np.linspace(1e-9,1.8,15)

# <codecell>
# Load the data
ne_sims = []
for ii in range(len(all_sims)):
    pluto_dir = os.path.join(os.path.expandvars(all_sims[ii]),'out')
    q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                # time=125.0,
                                                nframe=pluto_nframes[ii]
                                                )
    # Convert r and z to cm
    r /= 1e3
    z /= 1e3

    ne_integ_z = []
    for jj in range(len(z_lines)):
        # Integration on line z=const
        # z_lines will be between z[idx_z] and z[idx_z+1]
        idx_z = np.argmax(z[z_lines[jj]>=z])
        # Just to debug
        # print("z[idx_z]={};\tz_line={};\tz[idx_z+1]={}".format(z[idx_z],
        #                                                      z_lines[jj],
        #                                                      z[idx_z+1]))
        integ = np.sum(np.pi * q["ne"][idx_z,:] * (r[1:]**2 - r[:-1]**2))
        ne_integ_z.append(integ)

    ne_sims.append(np.array(ne_integ_z))

# <codecell>
# Plots
fig, ax = plt.subplots()

for ii in range(len(all_sims)):
    ax.plot(z_lines, ne_sims[ii], label=legends[ii])
ax.legend()
