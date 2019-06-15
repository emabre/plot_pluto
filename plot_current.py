import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import os

# <codecell>
# Settings
# sim = '$TORRE/home/konrad/simulazioni/sims_pluto/I90/newtransp'
# sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho8e-8-I235-3.2cmL-1mmD'
sim = '/home/ema/simulazioni/sims_pluto/perTesi/rho6e-7-I245-1.2cmL-1mmD-NEWGRID/'

# <codecell>
# Load data
current_finame = 'current_table.dat'
curr_fi = os.path.join(os.path.expandvars(sim),current_finame)

curr = np.loadtxt(curr_fi)

# <codecell>
# Plot
fig, ax = plt.subplots()
ax.plot(curr[:,0], curr[:,1], 'g')
ax.set_xlabel('Time / s')
ax.set_ylabel('Current / A')
ax.set_title(curr_fi)
ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
ax.grid()
plt.show()
