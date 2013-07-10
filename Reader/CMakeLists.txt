cmake_minimum_required(VERSION 2.8)

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules;${CMAKE_MODULE_PATH}"
)
FIND_PACKAGE(HDF5 REQUIRED)
FIND_PACKAGE(Wcs REQUIRED)
FIND_PACKAGE(CasaCore REQUIRED
 casa_fits
 casa_measures
 casa_scimath
 casa_components
 casa_images
 casa_mirlib
 casa_scimath_f
 casa_coordinates
 casa_lattices
 casa_ms
 casa_tables
 casa_derivedmscal
 casa_meas
 casa_msfits
 casa_casa
)
FIND_PACKAGE(ParaView REQUIRED)
INCLUDE(${PARAVIEW_USE_FILE})
message("VTK_REQUIRED_CXX_FLAGS: ${VTK_REQUIRED_CXX_FLAGS}")

include_directories("src")


IF(PARAVIEW_BUILD_QT_GUI)
  include_directories("${WCS_INCLUDE_DIR}")
  include_directories("${CASACORE_INCLUDE_DIR}")
  include_directories("${HDF5_INCLUDE_DIRS}")

  # The pqReader.xml file contains xml defining readers with their file
  # extensions and descriptions.  The Server Manager XML file defines the
  # properties of the reader.  This example uses the already existing PNG
  # reader but gives it another identity in the Server Manager.  Adding your
  # own reader to ParaView would mean the SourceProxy class points to your
  # class.

  ADD_PARAVIEW_PLUGIN(FitsReader "1.0"
    SERVER_MANAGER_XML server_manager.xml
    GUI_RESOURCE_FILES gui_resource.xml
    SERVER_MANAGER_SOURCES
      src/FitsReader.h src/FitsReader.cpp
      src/SliceImage.h src/SliceImage.cpp
      src/AccumulateSlices.h src/AccumulateSlices.cpp
      src/LofarGaussianBlur.h src/LofarGaussianBlur.cpp
      src/LofarMaskNoise.h src/LofarMaskNoise.cpp
      src/LofarApplyMask.h src/LofarApplyMask.cpp
      src/LofarIntegrateFrequencies.h src/LofarIntegrateFrequencies.cpp
    )
  target_link_libraries(FitsReader ${HDF5_LIBRARY} ${CASACORE_LIBRARIES} ${WCS_LIBRARIES})

#  add_executable(test_application src/main.cpp)
#  target_link_libraries(test_application FitsReader)

#  add_executable(regrid src/regrid.cpp)
#  target_link_libraries(regrid ${HDF5_LIBRARY} ${CASACORE_LIBRARIES} ${WCS_LIBRARIES})
ENDIF(PARAVIEW_BUILD_QT_GUI)

