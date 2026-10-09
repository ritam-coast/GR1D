F90=gfortran -std=legacy -fallow-argument-mismatch -ffree-line-length-none -ffixed-line-length-none
F90FLAGS= -g -O3 -fopenmp 
LDFLAGS= -g -O3 -fopenmp 
MODINC="-I ./"

HDF5DIR=/usr
HDF5INCS=-I/usr/include/hdf5/serial
HDF5LIBS=-L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5_fortran -lhdf5 -lz

LAPACKDIR=
LAPACKLIBS=-L${LAPACKDIR}/lib -llapack -lblas

HAVE_NUC_EOS=0
HAVE_RESTART=1
HAVE_LEAK_ROS=0
HAVE_LAPACK=0
