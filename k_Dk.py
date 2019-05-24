import numpy as np
import importlib
import os
import pluto_read_frm as prf
import utilities as ut

importlib.reload(prf)
importlib.reload(ut)

def k_Dk_time_evol(sim, pluto_nframes, r_cap, l_cap):
    '''
    pluto_nframes : the frames of pluto which I want to see (it must be a list of integers)
    r_cap : Capillary radius
    l_cap : Capillary length
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

    # Build the focusing strenght k(r)=B/r (vertex centered)
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

    # Build k (it is not really k, just B/R)
    k = np.zeros(B_avg_z_v_c.shape)
    for ii in range(B_avg_z_v_c.shape[1]):
        k[1:,ii] = B_avg_z_v_c[1:,ii]/r_c[1:]
    k[0,:] = k[1,:]
    # Build delta k
    Dk = np.zeros(k.shape)
    for ii in range(k.shape[1]):
        Dk[:,ii] = k[0,ii] - k[:,ii]

    return times, r_c, k, Dk
