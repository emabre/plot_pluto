import os
import numpy as np
def get_currtab(sim_path, current_finame = 'current_table.dat'):
    '''
    Returns time and current (as two separate 1D vectors of same size) of the current table
    '''
    curr_fi = os.path.join(os.path.expandvars(sim_path),current_finame)
    curr = np.loadtxt(curr_fi)
    return curr[:,0], curr[:,1]  # curr[:,0] is the time, curr[:,1] is the current.
