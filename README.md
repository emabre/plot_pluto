# plot_pluto
Some python and R scripts to plot the results of capillary discharge simulations made with the MHD code PLUTO.

## Prerequisites
Besides quite usual python modules (os, numpy, matplotlib..) you need **vtk**

## How to use it
Simply import pluto_read_frm module, than use the function of that module:
``q, r, z, theta, t, n = pluto_read_vtk_frame(path_to_vtk_file, nframe=frame_number)``
Inside q you have a dictionary, and its keys are the variable names (e.g. "bx3"), its values are numpy arrays.
