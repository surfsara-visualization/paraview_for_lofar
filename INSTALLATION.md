Installation
==================

The installation uses CMake and depends on some external libraries:
* paraview (including the source, I used 3.98.1)
  * Qt (at least version 4.7)
* casacore
  * wcs
  * hdf5
  * cfitsio3
  * lapack

Installation of paraview
==================
CMAKE_BUILD_TYPE: Release
PARAVIEW_INSTALL_DEVELOPMENT_FILES: ON
PARAVIEW_USE_MPI: ON
PARAVIEW_ENABLE_PYTHON: ON