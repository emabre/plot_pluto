import numpy as np
import matplotlib.pyplot as plt
import importlib
import pluto_read_frm as prf

importlib.reload(prf)
plt.close("all")

# <codecell>
pluto_dir = '/home/ema/simulazioni/sims_pluto/I90/newtransp-rho20'
pluto_nframe = 25
quantities = ['bx3']
ordering = 'F'

# <codecell>
# Load data
q, x, y, z, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                            # time=125.0,
                                            nframe=pluto_nframe
                                            )

# <codecell>
# Plotting
fig, ax = plt.subplots()

ax.pcolormesh(x, y, q["bx3"])

plt.show()
