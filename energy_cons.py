import numpy as np
import matplotlib.pyplot as plt

# <codecell>
en_cons_fi = "/home/ema/simulazioni/sims_pluto/disch_outcap/out/energy_cons.dat"

# <codecell>
en_cons = np.loadtxt(en_cons_fi)

t = en_cons[:, 0]
E_tc_in = en_cons[:, 5]

# <codecell>

fig, ax = plt.subplot()
ax.plot(t,en_cons)
