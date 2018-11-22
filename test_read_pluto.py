import numpy as np
import matplotlib.pyplot as plt
import importlib
import os
import capillary_utils as cu
import pluto_read_frm as prf

importlib.reload(prf)
importlib.reload(cu)

plt.close("all")

# <codecell>
sim = '/home/ema/simulazioni/sims_pluto/I90/newtransp-rho20'
pluto_nframe = 25
quantities = ['bx3']
ordering = 'F'

# <codecell>
# Load data
pluto_dir = os.path.join(sim,'out')
q, x, y, z, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                            # time=125.0,
                                            nframe=pluto_nframe
                                            )

# x,y = np.meshgrid(x,y)

# <codecell>
# Plotting
fig, ax = plt.subplots()

ax.pcolormesh(x, y, q["bx3"])

plt.show()
