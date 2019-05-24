import numpy as np
import importlib
import os
import pluto_read_frm as prf
import utilities as ut


importlib.reload(prf)
importlib.reload(ut)

def emittance(x, xp):
    cov_xxp = np.cov(np.stack([x,xp]))[0,1]
    sigma_x = x.std()
    sigma_xp = xp.std()
    emitt = np.sqrt(sigma_x**2 * sigma_xp**2 -  cov_xxp**2)
    return emitt, sigma_x, sigma_xp, cov_xxp

#
#def emittance_growth_thin_lens(x, xp, k):
#    
#    # Generate particle distribution

def g_Dg_time_evol(sim, pluto_nframes, r_cap, l_cap):
    '''
    Compute <B>/r and variation of <B>/r with respect to d<B>/dr(r=0).
    Input:
        pluto_nframes : the frames of pluto which I want to see (it must be a list of integers)
        r_cap : Capillary radius
        l_cap : Capillary length
    Returns:
        times, r_c, g, Dg
    '''

    # Load the data
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
        B.append(q["bx3"]/1e4)  # I convert B to Tesla

        # Build capillary shape (False where there is wall, True elsewere)
        # u_cap = cap_shape(r, z, r_cap, l_cap)
        cap.append(q['interBound'] == 0e0)

        # Averaging of B over z (I cannot use trapz, since I have cell-averaged B)
    #    B_integ_z = np.sum(((B[ii]*cap[ii].astype(int)).T*np.diff(z)).T,
    #                                       axis=0)
        B_integ_z = np.sum(((B[ii]*cap[ii]).T*np.diff(z)).T,
                                           axis=0)
        B_avg_z.append(2*B_integ_z/l_cap)

    # Build the focusing strenght g(r)=B/r (vertex centered)
    times = np.array(times)
    # Vertex centered integral of B
    B_avg_z_v = []
    for ii in range(len(B_avg_z)):
        B_avg_z_v.append(np.concatenate((np.array([0]),
                                         0.5*(B_avg_z[ii][1:] + B_avg_z[ii][:-1]),
                                         np.array([np.nan]))))

    # I radially restrict the B field to the capillary
    B_avg_z_v_c = np.transpose(np.array([list(B_avg_z_v[ii][r<=r_cap]) for ii in range(len(B_avg_z_v))]))
    r_c = r[r<=r_cap]

    # Build g (it is not really g, just B/R)
    g = np.zeros(B_avg_z_v_c.shape)
    for ii in range(B_avg_z_v_c.shape[1]):
        g[1:,ii] = B_avg_z_v_c[1:,ii]/r_c[1:]
    g[0,:] = g[1,:]
    # Build delta g
    Dg = np.zeros(g.shape)
    for ii in range(g.shape[1]):
        Dg[:,ii] = g[0,ii] - g[:,ii]

    return times, r_c, g, Dg
