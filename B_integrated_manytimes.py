import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import pluto_read_frm as prf

#plt.close("all")
importlib.reload(prf)
# plt.ion()

# <codecell>
# Settings
paper_emulate = 'Pompili2017'

if paper_emulate=='Pompili2017':
    pluto_nframes = [130]
    # Capillary radius
    r_cap = 0.5e-3
    # Capillary length, half of the real one
    l_cap = 1.5e-2
    sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho4.5e-7-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
elif paper_emulate==None:
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1e5Pa'
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa'
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-rhounif-I90-3.2cm'
    sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
    B_measured = np.loadtxt('/home/ema/Dottorato/dati_sperimentali_e_calcoli/Tabulazione_esperimAPL/ArticoloPompili2017/extracted_data/B.dat')
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho4.5e-7-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
    # Capillary length, half of the real one
    l_cap = 1.5e-3
    # Capillary radius
    r_cap = 0.5e-3
    pluto_nframes = [100,130]

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
    ax_avg.plot(0.5*(r[:-1]+r[1:])*1e6, B_avg_z[ii],
                '.-',
                label='simulated',
                )
# ax_avg.plot(B_measured[:,0], ,'.-')

ax_avg.legend()
ax_avg.set_title('t = {} ns'.format(times[ii]))
ax_avg.set_xlim([0.,r_cap*1e6])
ax_avg.set_ylabel('Mag. field, longitudinal average (T)')
ax_avg.set_xlabel('r (μm)')
ax_avg.grid()

plt.show()
#ax_zt.set_aspect(1)
