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
paper_emulate = 'Pompili2017'

if paper_emulate=='Pompili2017':
    pluto_nframes = [130]
    # Capillary radius
    r_cap = 0.5e-3
    # Capillary length, half of the real one
    l_cap = 1.5e-2
    # Measured B as in Pompili2017
    B_measured = np.loadtxt('/home/ema/Dottorato/dati_sperimentali_e_calcoli/Tabulazione_esperimAPL/ArticoloPompili2017/extracted_data/B.dat')
    # Simulated with pluto
    sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho4.5e-7-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
    # Static B as in Bobrova model
    B_staticBob = np.loadtxt('/home/ema/myprogr/scarica_a_regime/B_adimensional.dat')
elif paper_emulate==None:
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1e5Pa'
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa'
    # sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-rhounif-I90-3.2cm'
    sim = '/home/ema/simulazioni/sims_pluto/dens_real/1.3e5Pa-1.2cm'
    # sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho4.5e-7-I90-3.2cmL-1mmD-r60-NTOT16-diffRecPeriod8'
    # Capillary length, half of the real one
    l_cap = 1.5e-2
    # Capillary radius
    r_cap = 0.5e-3
    pluto_nframes = [100,130]

# <codecell> Load the data
B = []
r_cc = []; z_cc = []
cap = []
times = []
B_avg_z = []; ioniz_z0 = []
for ii in range(len(pluto_nframes)):
    pluto_dir = os.path.join(os.path.expandvars(sim), 'out')
    q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                    # time=125.0,
                                                    nframe=pluto_nframes[ii])
    times.append(t)
    # Convert r and z to m
    r /= 1e5
    z /= 1e5
    r_cc.append(0.5*(r[:-1]+r[1:]))
    z_cc.append(0.5*(z[:-1]+z[1:]))

    B.append(q["bx3"]*1e-4)  # Convert to Tesla

    # Build capillary shape (False where there is wall, True elsewere)
    # u_cap = cap_shape(r, z, r_cap, l_cap)
    cap.append(q['interBound'] == 0e0)

    # Averaging of B over z (I cannot use trapz, since I have cell-averaged B)
    B_integ_z = np.sum(((B[ii]*cap[ii].astype(int)).T*np.diff(z)).T, axis=0)
    B_avg_z.append(B_integ_z/(l_cap))

    ioniz_z0.append(q["ioniz"][np.argmin(np.abs(z_cc[ii])),:])

# <codecell>
# Plots
# Average ne on fixed z positions
fig_avg, ax_avg = plt.subplots(figsize=(3.5, 2.7))
# Fix artifically the last point of B, otherwise bad looking (due to mesh)
B_avg_z[0][np.argmin(np.abs(500.-0.5*(r[:-1]+r[1:])*1e6))+1] = 1e-3*B_measured[-1,1]
# plot B computed and measured
ax_ioniz = ax_avg.twinx()
ax_ioniz.set_zorder(1); ax_avg.set_zorder(10)
ax_avg.patch.set_visible(False)
ax_ioniz.tick_params(axis='y', labelcolor='olivedrab')
for ii in range(len(pluto_nframes)):
    ax_ioniz.plot(r_cc[ii]*1e6, ioniz_z0[ii],
                '--', color='olivedrab',
                label='ioniz. degree',
                # zorder = 1,
                )
    ax_avg.plot(r_cc[ii]*1e6, B_avg_z[ii]*1e3,
                '-', color='purple',
                label='Simulated',
                # zorder = 10,
                )
ax_ioniz.set_ylim(0.,1.)
ax_ioniz.set_ylabel('Ionization degree')
ax_avg.plot(B_measured[:,0], B_measured[:,1],
            '-',
            color='k',
            label='Exp. inferred')
ax_avg.plot(B_staticBob[:,0]*r_cap*1e6/B_staticBob[-1,0],
            B_staticBob[:,1]*B_measured[-1,1]/B_staticBob[-1,1],
            '-',
            color='darkorange',
            label='Equil. model')

ax_avg.legend(framealpha = 0.4).set_zorder(102)
ax_avg.set_title('t = {:.0f}ns, I = {:g}A'.format(times[ii],
                                            np.interp(times[ii]*1e-9, *ut.get_currtab(sim))))
ax_avg.set_xlim([0.,r_cap*1e6])
ax_avg.set_ylim([0., 25.])
ax_avg.set_ylabel('Magnetic field (longitud. avg.) (mT)')
ax_avg.set_xlabel('r (μm)')
ax_avg.grid()
fig_avg.tight_layout()
plt.show()
#ax_zt.set_aspect(1)
