#!/usr/bin/python3
# -*-coding:utf-8 -*


from mpi4py import MPI
mpi_comm = MPI.COMM_WORLD
mpi_rank = mpi_comm.Get_rank()
mpi_size = mpi_comm.Get_size()
mpi_name = mpi_comm.Get_name()

from rayleigh_utils import GlobalTimer

import cython
import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np

# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.double

# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.double_t DTYPE_t


# @GlobalTimer
def FastInterp(double x, np.ndarray[DTYPE_t] x_list, np.ndarray[DTYPE_t] y_list):
    cdef double interp = np.interp(x, x_list, y_list)
    return interp
