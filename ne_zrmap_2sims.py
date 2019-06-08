import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import importlib
import os
import pluto_read_frm as prf
from scipy.constants import mu_0
import utilities as ut
import active_plasma_lens as apl
import matplotlib.image as mpimg
from scipy.interpolate import griddata

#plt.close("all")
importlib.reload(prf)
importlib.reload(ut)
importlib.reload(apl)

# <codecell>
# Settings

# ----
sim = ['/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I500flattop-1.2cmL-1mmD',
       '/home/ema/simulazioni/sims_pluto/perTesi/rho2.53e-7-I720flattop-1.2cmL-1.2mmD'
       ]
pluto_nframe = 50
l_cap = 1e-2 # m
r_cap = 0.5e-2

#%% Load the data
if len(sim)!=2:
    raise ValueError('sim must be a list of two(2) simulation paths')
ne = [[], []]
r = [[],[]]
z = [[],[]]
t = [[],[]]
for ss in range(len(sim)):
    pluto_dir = os.path.join(os.path.expandvars(sim[ss]), 'out')
    q, r[ss], z[ss], theta, t[ss], n = prf.pluto_read_vtk_frame(pluto_dir,
                                                                # time=125.0,
                                                                nframe=pluto_nframe)
    ne[ss] = q["ne"]
    # Convert r and z to m
    r[ss] /= 1e5
    z[ss] /= 1e5

#%%
# Interpolate with nearest the border points
r_cc = [0.5*(r[ss][1:]+r[ss][:-1]) for ss in (0,1)]
z_cc = [0.5*(z[ss][1:]+z[ss][:-1]) for ss in (0,1)]

rr_cc = [[],[]]; zz_cc = [[],[]]
for ss in (0,1):
    rr_cc[ss], zz_cc[ss] = np.meshgrid(r_cc[ss], z_cc[ss])

points = [[],[]]; ne_cc = [[],[]]
for ss in (0,1):
    for ii in range(len(r_cc[ss])):
        for jj in range(len(z_cc[ss])):
            # points.append([rr_cc[ss], zz_cc[ss]])
            points[ss].append([z_cc[ss][jj], r_cc[ss][ii]])
            ne_cc[ss].append(ne[ss][jj,ii])
ne_cc = [np.array(ne_cc[ss]) for ss in range(len(ne_cc))]
points = [np.array(points[ss]) for ss in range(len(points))]


zz = [[],[]]; rr = [[],[]]
for ss in (0,1):
    zz[ss], rr[ss] = np.meshgrid(z[ss],r[ss])

# vertex centered ne
ne_vc = [griddata(points[0], ne_cc[0],
                  np.stack((zz[ss].flatten(), rr[ss].flatten()), axis=1),
                  method='linear') for ss in (0,1)]
ne_vc = [np.reshape(ne_vc[ss], (len(z[ss]),len(r[ss]))) for ss in (0,1)]


fig, ax = plt.subplots()
ss = 0
ax.contourf(r[0], z[0], ne_vc[0])
