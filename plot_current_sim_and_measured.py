import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import os


#%% Settings

# Allowed choices are '245A-1cm' and '90A-3cm'
measure_choice = '245A-1cm'  # '90A-3cm' or '245A-1cm'

# Imposed current
if measure_choice=='245A-1cm':
    sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho6e-7-I245-1.2cmL-1mmD-NEWGRID'
    # Measured current
    disch_data  = os.path.join(os.path.expandvars('$DOTTORATO'),
                                  'dati_sperimentali_e_calcoli',
                                  'misure_capillare_1cmL-1mmD',
                                  'stessa_scala_temporale',
                                  'corrente',
                                  'Current_ASCII.txt')
    Imax = 250.
elif measure_choice=='90A-3cm':
    # Measured current
    # sim = '$TORRE/home/konrad/simulazioni/sims_pluto/I90/newtransp'
    sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-6-I90-3.2cmL-1mmD'
    disch_data  = os.path.join(os.path.expandvars('$DOTTORATO'),
                                 'dati_sperimentali_e_calcoli',
                                 'misure_capillare_3cmL-1mmD',
                                 'stessa_scala_temporale',
                                 'corrente',
                                 'Current_ASCII.txt')
    Imax = 100.
# Only times higher than t_start_s (seconds) will be plotted of the measured data
t_start_s = 0.  # -30e-9
# Only times lower than t_end_s (seconds) will be ploted of the measured data
t_end_s = 2e-6
# Only one point every nskip points will be plotted
nskip = 10
# Settings of min max time for plot (s)
t_min = 0.
t_max = 2e-6

#%% Load data
# measured data
disch = np.loadtxt(disch_data, skiprows=1)
t = disch[:,0]
I = disch[:,1]
# simlation data
current_finame = 'current_table.dat'
curr_fi = os.path.join(os.path.expandvars(sim),current_finame)
curr = np.loadtxt(curr_fi)
t_sim = curr[:,0]
I_sim = curr[:,1]

#%% Plot
fig, ax = plt.subplots(figsize=(3.85,2.7))

ax.plot(t[(t>t_start_s) & (t<t_end_s)][::nskip]*1e9,
        I[(t>t_start_s) & (t<t_end_s)][::nskip], 'g',
        label='Measured')
ax.plot(t_sim*1e9, I_sim, '--', color='k',
        label='Set in sim.')
ax.set_ylim([0.,Imax])
ax.set_xlim([t_min*1e9, t_max*1e9])
ax.legend()

# ax[0].set_title(r'\SI{2}{\nano\farad}-\SI{100}{\ohm} circuit')
# ax[1].set_title(r'\SI{7.2}{\nano\farad}-\SI{36.5}{\ohm} circuit')
# ax.set_title('$2\mathrm{nF}-100\Omega$ circuit,\n$3\mathrm{cm}$ length capillary')
# ax.set_title('$7.2\mathrm{nF}-36.5\Omega$ circuit,\n$1\mathrm{cm}$ length capillary')

ax.set_ylabel('Current (A)')
ax.set_xlabel('Time (ns)')
ax.grid()
# ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))

fig.tight_layout()
plt.show()
