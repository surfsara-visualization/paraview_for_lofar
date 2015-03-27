Paraview for Lofar
==================

The repository contains a paraview reader for Lofar data in the FITS or HDF5 format, using the casacore library.

There are also several filters to manipulate the data:
* AccumulateSlices<br>
* LofarApplyMask<br>
  Given an image an a mask, set a constant value to the masked region
* LofarGaussianBlur<br>
  Apply a Gaussian blur to each of the frequency images separately
* LofarIntegrateFrequencies<br>
  Integrate the images along the frequency direction (z-axis)
* LofarMaskNoise<br>
  Create a mask for every frequency image where the signal is less than a given signal, given in terms of the standard deviation of the frequency image.
* SliceImage<br>
  Extracts a slice through the dataset that is orthogonal to one of the axis. The slice can be along a line segment or a spline.

Finally, it contains a Lofar-specific ParaView view called LofarView 