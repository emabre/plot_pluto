# import pluto_read_frm as prf
import numpy as np
import matplotlib.pyplot as plt
import importlib
import pluto_read_frm as prf
# from vtk import vtkRectilinearGridReader
# from vtk.util import numpy_support as VN
importlib.reload(prf)

# <codecell>
# pluto_vtk =  '/home/ema/simulazioni/sims_pluto/I90/newtransp-rho20/out/data.0025.vtk'
pluto_dir = '/home/ema/simulazioni/sims_pluto/I90/newtransp-rho20'
pluto_nframe = 25
quantities = ['bx3']
ordering = 'F'

# <codecell>
# Load data
q, x, y, z = prf.pluto_read_vtk_frame(pluto_dir, pluto_nframe)
# u, x, y, z = prf.pluto_vtk_to_numpy(pluto_vtk, quantities, ordering)

# <codecell>
# Plotting
fig, ax = plt.subplots()

ax.pcolormesh(x, y, q["bx3"])
