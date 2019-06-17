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

def generate_beam_transverse(sigma_x, d_sigma_x, emitt_x, Npart, check_distro=True):
    '''Build transverse distribution of particles inside beam from rms emittance,
    rms size and z-derivative of rms size.'''

    cov_xxp = sigma_x * d_sigma_x
    sigma_xp = np.sqrt((emitt_x**2 + cov_xxp**2)/sigma_x**2)
    # mean and cov
    mean = [0,0]
    cov = [[sigma_x**2, cov_xxp], [cov_xxp, sigma_xp**2]]
    # Generate distribution
    x, xp = np.random.multivariate_normal(mean, cov, Npart).T

    if check_distro:
        print('generated distro with {} partilces'.format(Npart))
        # print(emitt_x)
        print('Emittance relative error: {}'.format((emitt_x-emittance(x, xp)[0])/emitt_x))
        print('Spot:')
        print('Spot rms, required:{:.5g}; obtained: {:.5g}'.format(sigma_x,
                                                                   np.std(x)))
    return x, xp

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

def ne_avg_over_r(sim, pluto_nframes, average_ne, z_lines=None, ret_z_cell_borders=False):
    '''
    Compute max electron density over the transverse section (r direction).
    return a vector of ne, the elements of ne correspond to different z positions
    Input:
        pluto_nframes : the frames of pluto which I want to see (it must be a list of integers)
    Returns:
        times, e
    '''

    times = []
    ne_avg_r_sims = []
    for ii in range(len(pluto_nframes)):
        pluto_dir = os.path.join(os.path.expandvars(sim),'out')
        # Load the data
        q, r, z, theta, t, n = prf.pluto_read_vtk_frame(pluto_dir,
                                                        # time=125.0,
                                                        nframe=pluto_nframes[ii])
        times.append(t)
        # Convert r and z to cm
        r /= 1e3
        z /= 1e3
        ne_avg_r = []
        # If z_lines have not been provided, I take all z cell centers
        if z_lines is None:
            # idx_z_all = list(range(len(z)-1))  # I reduce by 1, since z are the cell borders, not the centers
            z_lines = 0.5*(z[1:]+z[:-1])
        for z_line in z_lines:
            # Check that z-lines are inside the grid
            if z_line>z.max() or z_line<z.min():
                raise ValueError('z_lines outside z grid')
            # I find the grid cell where z_line is included,
            # note that ne is defined inside the cell (it's an average value)
            # so I may imagine that it is constant inside the cell
            idx_z = np.argmax(z[z<=z_line])
            if average_ne == 'integral' or average_ne == 'mean' :
                # Build capillary shape (False where there is wall, True elsewere)
                cap = (q['interBound']==0e0)
                areas = []
                integ = np.sum(np.pi * q["ne"][idx_z,:] * (r[1:]**2 - r[:-1]**2) * cap[idx_z,:])
                area_r = np.sum(np.pi * (r[1:]**2 - r[:-1]**2) * cap[idx_z,:])
                areas.append(area_r)
                ne_avg_r.append(integ/area_r)
            elif average_ne == 'max':
                ne_avg_r.append(q["ne"][idx_z,:].max())
            else:
                raise ValueError('Wrong choice for average_ne')
        ne_avg_r_sims.append(ne_avg_r)
    if ret_z_cell_borders:
        return z_lines, z, ne_avg_r_sims, times
    else:
        return z_lines, ne_avg_r_sims, times
