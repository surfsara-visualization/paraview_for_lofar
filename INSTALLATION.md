Dependencies
============

The installation uses CMake and depends on some external libraries:
* ParaView (including the source, I used 3.98.1)
  * Qt (at least version 4.7)
* casacore
  * wcs
  * hdf5
  * cfitsio3
  * lapack

CMake flags to use when building paraview
=========================================

* CMAKE_BUILD_TYPE: Release
* PARAVIEW_INSTALL_DEVELOPMENT_FILES: ON
* PARAVIEW_USE_MPI: ON
* PARAVIEW_ENABLE_PYTHON: ON

Building the plugins
====================

1. Create a build directory, e.g build/
2. cd build/
3. $ ccmake -DParaView_DIR=<pv-3.98.1-dir>/lib/cmake/paraview-3.98 ../Code/
4. Set LOFAR_BUILD_FILTERS and LOFAR_BUILD_LOFARVIEW as appropriate
5. Set other CMake options as appropriate
6. [C]onfigure, [G]enerate
7. make


LofarView
=========

The LofarView plugin will be called libLofarView.so and will be located
in the build directory.

Installing the plugin in ParaView
---------------------------------

1. Start ParaView 3.98.1
2. Tools -> Manage Plugins
3. Under "Local Plugins" pick "Load New..."
4. Navigate to the directory where the plugin was built
5. Choose the file libLofarView.so
6. Press Ok
   The plugin list should now show LofarView as Loaded"
7. Close the plugin manager window

Creating a LofarView view
-------------------------

1. Close the main 3D View (with the X in the upper right)
2. Pick "LofarView" in the "Create View" list

If no data is loaded the lower-right of the view should show "No Data".
Whenever the mouse is clicked somewhere within the view the cursor
coordinates shown will update.
readme.txt (END)


