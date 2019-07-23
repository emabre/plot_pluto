import numpy as np
import importlib
import os
import pluto_read_frm as prf
import utilities as ut
import scipy.constants as cst


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

def g_Dg_time_evol(sim, pluto_nframes, r_cap, l_cap, ret_full_g=False):
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

        # Averaging of B over z (I cannot use trapz, since I have cell-averaged B)
        # Beware: also the bcs on B (i.e. the first cells near the capillar wall that are in the internal boundary)
        # are useful for computing, because if I want to compute B at the cell edges, I need them
        B_integ_z = np.sum((B[ii].T*np.diff(z)).T, axis=0)

        B_avg_z.append(2*B_integ_z/l_cap)

    # Build the focusing strenght g(r)=B/r (vertex centered)
    times = np.array(times)
    # Vertex centered averaged B
    B_avg_z_v = []
    B_v = []
    for ii in range(len(pluto_nframes)):
        B_avg_z_v.append(np.concatenate((np.array([0]),
                                         0.5*(B_avg_z[ii][1:] + B_avg_z[ii][:-1]),
                                         np.array([np.nan]))))
        B_v.append(np.concatenate((np.zeros((B[ii].shape[0]-1,1)),
                                   0.25*(B[ii][1:,1:]+B[ii][1:,:-1]+B[ii][-1:,1:]+B[ii][-1:,-1])),
                                  axis = 1))
    z_B_v = z[1:-1]
    r_B_v = r[:-1]

    # I radially restrict the B field to the capillary
    B_avg_z_v_c = np.transpose(np.array([list(B_avg_z_v[ii][r<=r_cap]) for ii in range(len(B_avg_z_v))]))
    # B_v = np.transpose(np.array([list(B_v[ii][r<=r_cap]) for ii in range(len(B_v))]))
    r_c = r[r<=r_cap]

    # Build the longitudinal average of g (it is not really g, just B/R)
    g = np.zeros(B_avg_z_v_c.shape)
    for ii in range(B_avg_z_v_c.shape[1]):
        g[1:,ii] = B_avg_z_v_c[1:,ii]/r_c[1:]
    g[0,:] = g[1,:]
    # Build delta g
    Dg = np.zeros(g.shape)
    for ii in range(g.shape[1]):
        Dg[:,ii] = g[0,ii] - g[:,ii]

    # Build g as 2D matrix (to represent a function of (r,z))
    g_full = []
    for ii in range(len(pluto_nframes)):
        g_temp = B_v[ii][:,1:]/r_B_v[1:]
        g_full.append(np.concatenate((g_temp[:,[0]],
                                      g_temp),
                                      axis=1))


    if ret_full_g:
        return times, r_c, g, Dg, B_v, g_full, z_B_v, r_B_v
    else:
        return times, r_c, g, Dg

def ne_avg_over_r(sim, pluto_nframes, average_ne, z_lines=None, ret_z_cell_borders=False):
    '''
    Compute max/integral-avg/axial electron density over the transverse section (r direction).
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
            elif average_ne == 'axis':
                ne_avg_r.append(q["ne"][idx_z,0])
            else:
                raise ValueError('Wrong choice for average_ne')
        ne_avg_r_sims.append(ne_avg_r)
    if ret_z_cell_borders:
        return z_lines, z, ne_avg_r_sims, times
    else:
        return z_lines, ne_avg_r_sims, times


def focus_in_thin_apl(g, r_c, x, xp, y, l_cap, gamma, Dz):
    '''
    Focus a beam passing through an APL in thin lens approximation.
    g: mag field gradient (B/R, Tesla/m) (1D array like)
    r_c: radial points where g is defined (m) (1D array like)
    x: transverse beam particle positions (m) (1D array like)
    xp: angular divergence of beam particles (m) (1D array like)
    l_cap: capillary length (m)
    gamma: beam relativistic gamma
    Dz: drift after lens (m)
    Returns
    sigma_x_new,
    emitt_x_new, (non normalized emittance after lens)
    x_new,
    xp_new
    '''
    if len(x)!=len(y):
        raise ValueError('x and y must have same length')
    if len(r_c)!=len(g):
        raise ValueError('r_c and g must have same length')

    # Interpolate field gradient at the particle positions
    g_real_interp = np.interp(np.sqrt(x**2+y**2),
                              np.concatenate((np.flip(-r_c[1:], axis=0), r_c)),
                              np.concatenate((np.flip(g[1:], axis=0), g)))

    # Define focusing strength
    k = cst.e/(cst.m_e*cst.c*gamma) * g_real_interp
    K = k*l_cap

    # Divergence increase in thin lens approx
    Dxp = np.zeros(len(K))
    xp_new = np.zeros(len(K))
    Dxp = - K*x
    xp_new = xp + Dxp

    # New emittance after lens
    emitt_x_new = emittance(x, xp_new)[0]

    # New spot after drift Dz following lens
    x_new = x + Dz*(xp_new)
    sigma_x_new = np.std(x_new)

    return sigma_x_new, emitt_x_new, x_new, xp_new

def focus_in_thick_apl_new(g, r_c, x, xp, y, yp, l_cap, gamma, Dz, Nz = 100):
    '''
    Focus a beam passing through an APL as thick lens. The beam is assumed to have zero thickness in z direction (longitudinal),
    this has no effect on the tracking since the space charge and the wakefields are neglected (and we are not interested in the
    beam properties inside the lens but only outside).
    g: mag field gradient (B/R, Tesla/m) (1D array like)
    r_c: radial points where g is defined (m) (1D array like)
    x: transverse beam particle positions (m) (1D array like)
    xp: angular divergence of beam particles (m) (1D array like)
    l_cap: capillary length (m)
    gamma: beam relativistic gamma
    Dz: drift after lens (m)
    Returns
    sigma_x_new,
    emitt_x_new, (non normalized emittance after lens)
    x_new,
    xp_new
    '''
    if len(x)!=len(y):
        raise ValueError('x and y must have same length')
    if len(r_c)!=len(g):
        raise ValueError('r_c and g must have same length')

    from scipy.interpolate import griddata

    # Reflect g and r across r=0
    r_c_refl = np.concatenate((np.flip(-r_c[1:], axis=0), r_c))
    no
    g_refl = np.concatenate((np.flip(g[1:,:], axis=0), g))

    # I do a leapfrog
    z = 0; dz = l_cap/Nz
    # x_new_out = []; xp_new_out = [];
    # y_new_out = []; yp_new_out = []
    x_old = np.copy(x); y_old = np.copy(y)
    xp_old = np.copy(xp); yp_old = np.copy(yp)
    ii =0
    while z<l_cap:
        ii += 1
        # print('step {}'.format(ii))
        # Interpolate field gradient at the new particle positions
        # g_real_interp = np.interp(np.sqrt(x_old**2+y_old**2),
        #                           r_c_refl,
        #                           g_refl)
        # rz_beam = qualcosa da: np.sqrt(x_old**2+y_old**2), z
        g_real_interp = griddata(rz_c_refl, g_refl, rz_beam, method='linear')

        # Define focusing strength experienced by each particle (g_real_interp has been interpolated at particle positions)
        k = cst.e/(cst.m_e*cst.c*gamma) * g_real_interp

        # Divergence (xp,yp) increse
        xp_new = xp_old - k*dz*x_old
        yp_new = yp_old - k*dz*y_old

        # Update positions
        x_new = x_old + dz*xp_new
        y_new = y_old + dz*yp_new

        # Backup position and divergences for next iteration
        x_old = np.copy(x_new)
        y_old = np.copy(y_new)
        xp_old = np.copy(xp_new)
        yp_old = np.copy(yp_new)

        z += dz

    # # Save output data
    # x_new_out.append(x_new)
    # xp_new_out.append(xp_new)
    # y_new_out.append(y_new)
    # yp_new_out.append(yp_new)

    # New emittance after lens
    emitt_x_new = emittance(x_new, xp_new)[0]

    # New spot after drift Dz following lens
    x_new_drift = x_new + Dz*(xp_new)
    sigma_x_new = np.std(x_new_drift)
    y_new_drift = y_new + Dz*(yp_new)
    sigma_y_new = np.std(y_new_drift)

    return sigma_x_new, emitt_x_new, x_new_drift, xp_new, y_new_drift, yp_new

def focus_in_thick_apl(g, r_c, x, xp, y, yp, l_cap, gamma, Dz, Nz = 100):
    '''
    Focus a beam passing through an APL as thick lens.
    g: mag field gradient (B/R, Tesla/m) (1D array like)
    r_c: radial points where g is defined (m) (1D array like)
    x: transverse beam particle positions (m) (1D array like)
    xp: angular divergence of beam particles (m) (1D array like)
    l_cap: capillary length (m)
    gamma: beam relativistic gamma
    Dz: drift after lens (m)
    Returns
    sigma_x_new,
    emitt_x_new, (non normalized emittance after lens)
    x_new,
    xp_new
    '''
    if len(x)!=len(y):
        raise ValueError('x and y must have same length')
    if len(r_c)!=len(g):
        raise ValueError('r_c and g must have same length')

    # Reflect g and r across r=0
    r_c_refl = np.concatenate((np.flip(-r_c[1:], axis=0), r_c))
    g_refl = np.concatenate((np.flip(g[1:], axis=0), g))

    # I do a leapfrog
    z = 0; dz = l_cap/Nz
    # x_new_out = []; xp_new_out = [];
    # y_new_out = []; yp_new_out = []
    x_old = np.copy(x); y_old = np.copy(y)
    xp_old = np.copy(xp); yp_old = np.copy(yp)
    ii =0
    while z<l_cap:
        ii += 1
        # print('step {}'.format(ii))
        # Interpolate field gradient at the new particle positions
        g_real_interp = np.interp(np.sqrt(x_old**2+y_old**2),
                                  r_c_refl,
                                  g_refl)
        # Define focusing strength experienced by each particle (g_real_interp has been interpolated at particle positions)
        k = cst.e/(cst.m_e*cst.c*gamma) * g_real_interp

        # Divergence (xp,yp) increse
        xp_new = xp_old - k*dz*x_old
        yp_new = yp_old - k*dz*y_old

        # Update positions
        x_new = x_old + dz*xp_new
        y_new = y_old + dz*yp_new

        # Backup position and divergences for next iteration
        x_old = np.copy(x_new)
        y_old = np.copy(y_new)
        xp_old = np.copy(xp_new)
        yp_old = np.copy(yp_new)

        z += dz

    # # Save output data
    # x_new_out.append(x_new)
    # xp_new_out.append(xp_new)
    # y_new_out.append(y_new)
    # yp_new_out.append(yp_new)

    # New emittance after lens
    emitt_x_new = emittance(x_new, xp_new)[0]

    # New spot after drift Dz following lens
    x_new_drift = x_new + Dz*(xp_new)
    sigma_x_new = np.std(x_new_drift)
    y_new_drift = y_new + Dz*(yp_new)
    sigma_y_new = np.std(y_new_drift)

    return sigma_x_new, emitt_x_new, x_new_drift, xp_new, y_new_drift, yp_new

def drift_beam(x, xp, z_distances):
    '''
    Propagate a beam through a drift, and compute sigma and positions ad certain points,
    defined by the elements inside 1-D, array like Dz
    '''
    x_new = tuple(map(lambda dz: x+dz*xp, z_distances))
    sigma_x_new = np.array(tuple(map(np.std, x_new)))
    return x_new, sigma_x_new
