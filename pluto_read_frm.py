import numpy as np
# from vtk import vtkStructuredPointsReader
# from vtk import vtkStructuredGridReader
from vtk import vtkRectilinearGridReader
from vtk.util import numpy_support as VN
import os


def pluto_vtk_to_numpy(filename, quantity_names, ordering):
    '''pluto_vtk_to_numpy(filename, quantity_names, ordering): \
    loads a vtk as built by pluto \
    (only for scalar output of quantities)\
    and converts it to numpy array.\
    "ordering may be either "C" or "F"'''
    reader = vtkRectilinearGridReader()
    reader.SetFileName(filename)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()

    data = reader.GetOutput()

    dim = []
    for ii in data.GetDimensions():
        if ii > 1:
            dim.append(ii)

    # print('dim={}'.format(dim))

    # I transform to numpy array the domain (and reshape)
    x = np.zeros(data.GetNumberOfPoints())
    y = np.zeros(data.GetNumberOfPoints())
    z = np.zeros(data.GetNumberOfPoints())

    # Even though the data in pluto's vtk are defined on cells rather than points,
    # it is useful to have the points (which are the corners of the cells)
    for ii in range(data.GetNumberOfPoints()):
            x[ii], y[ii], z[ii] = data.GetPoint(ii)

    # print('x points: {}'.format(len(x)))
    # print('y points: {}'.format(len(y)))
    # print('z points: {}'.format(len(z)))

    x = x.reshape(dim, order=ordering).transpose()
    y = y.reshape(dim, order=ordering).transpose()
    z = z.reshape(dim, order=ordering).transpose()

    u = []
    for qn in quantity_names:
        # I transform to numpy array the scalar/vector field (and reshape)
        u.append(VN.vtk_to_numpy(data.GetCellData().GetArray(qn)))
        # # f = VN.vtk_to_numpy(data.GetField())
        # VN.vtk_to_numpy(data.GetPointData().GetArray(quantity_name))
    # field = data.GetFieldData()
    # time = VN.vtk_to_numpy(field.GetAbstractArray("time"))
    # current = VN.vtk_to_numpy(field.GetAbstractArray("current"))

    # Number of cells per each dimension
    dim_cells = list(np.array(dim)-1)
    dim_quantity = dim_cells.copy()
    for ii in range(len(u)):
        try:
            dim_quantity.append(u[ii].shape[1])
        except IndexError:
            pass
        u[ii] = u[ii].reshape(dim_quantity, order=ordering).transpose()

    return u, x, y, z

def pluto_read_vtk_frame(pluto_dir, n=None, t=None, q_names=None):
    '''Reads a vtk dataframe which was written by pluto 4.2, only in non-vector mode.
    Returns a dictionary with the quantities, ND arrays of x,y,z positions and an int or a float
    containig the actual time and number of frame'''

    # if n!=None and t==None:
    #     pluto_nframe = n
    #     pluto_time = leggi il tempo del frame n
    # elif n==None and t!=None:
    #     trovare l'n che più si avvicina a t
    #     pluto_nframe = qualcosa
    #     pluto_time
    # else
    #     raise ValueError("Epecify either n or t, not both or none.")

    # Names of the quantities (if not feeded as input)
    if q_names==None:
        q_names = ("rho","prs","vx1","vx2","ne","knor","ioniz","interBound",
                   "etax1", "c2p_fail", "bx3", "T", "Jz", "Jr")

    # Path of the file to read
    vtk_basename = "data."
    ordering = "F"
    vtk_finame = vtk_basename + '{:04d}.vtk'.format(pluto_nframe)
    vtk_fi = os.path.join(pluto_dir,"out",vtk_finame)

    # Read
    u, x, y, z = pluto_vtk_to_numpy(vtk_fi, q_names, ordering)
    # Build the dictionary
    q = {q_names[ii]:u[ii] for ii in range(len(q_names))}

    return q, x, y, z, pluto_time, pluto_nframe

def read_vtk_log(log_fi):
    '''Function to read pluto's vtk log files (vtk.out)'''

    with open(log_fi) as log:
        lines = log.readlines()

    nvtk = []; t = []; dt = []; nsteps = []
    file_type = []; endianess = []
    quantity_names = []
    for line in lines:
        elements = line.strip().split()
        nvtk.append(int(elements[0]))
        t.append(float(elements[1]))
        dt.append(float(elements[2]))
        nsteps.append(int(elements[3]))
        file_type.append(elements[4])
        endianess.append(elements[5])
        quantity_names.append(elements[6:])

    return (np.array(nvtk),
            np.array(t),
            np.array(dt),
            np.array(nsteps),
            file_type,
            endianess,
            quantity_names)

if __name__ == "__main__":
    pluto_vtk =  '/home/ema/simulazioni/sims_pluto/I90/newtransp-rho20/out/data.0025.vtk'
    quantities = ['bx3']
    u, x, y, z = pluto_vtk_to_numpy(pluto_vtk, quantities, 'C')
    len(u)
    x.shape
    y.shape
    u[0].shape
