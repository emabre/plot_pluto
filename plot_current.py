import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import os

# <codecell>
# Settings
# sim = '$TORRE/home/konrad/simulazioni/sims_pluto/I90/newtransp'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-8-I235-3.2cmL-1mmD'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho6e-7-I245-1.2cmL-1mmD-NEWGRID/'
sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD-NEWGRID/',
       '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1mmD-NEWGRID/']
legend = ['(a)','(b)']

# <codecell>
# Load data
current_finame = 'current_table.dat'
curr_fi = []; curr = []
for s in sim:
    curr_fi.append(os.path.join(os.path.expandvars(s),current_finame))
    curr.append(np.loadtxt(curr_fi[-1]))

# <codecell>
# Plot
fig, ax = plt.subplots(figsize=(3.85,2.7))
colors = ['g','crimson']
for ss in range(len(sim)):
    ax.plot(curr[ss][:,0]*1e9, curr[ss][:,1],
            c=colors[ss], label=legend[ss])
ax.set_xlabel('Time (ns)')
ax.set_ylabel('Current (A)')
ax.set_xlim(0,500)
ax.set_ylim(0,800)
# ax.set_title(curr_fi)
# ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
ax.legend()
ax.grid()
fig.tight_layout()
plt.show()
