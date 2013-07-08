# Locate casa-core
# This module defines
# CASACORE_LIBRARIES
# CASACORE_FOUND, if false, do not try to link to libcasacore 
# CASACORE_INCLUDE_DIR, where to find the headers
#
# $CASACORE_DIR is an environment variable that would
# correspond to the ./configure --prefix=$CASACORE_DIR
# used in building libcasacore.

set(_casacore_modules_to_process)
foreach(_cc_component ${CasaCore_FIND_COMPONENTS})
    list(APPEND _casacore_modules_to_process ${_cc_component})
endforeach()
list(APPEND _casacore_modules_to_process "dl" "cfitsio" "lapack")
list(REMOVE_DUPLICATES _casacore_modules_to_process)

FIND_PATH(CASACORE_INCLUDE_DIR fits/FITS.h
    NO_DEFAULT_PATH
    PATHS
    $ENV{CASACORE_DIR}/include
    $ENV{CASACORE_DIR}/include/casacore
    $ENV{CASACORE_DIR}
    ~/Library/Frameworks
    /Library/Frameworks
    /usr/local/include
    /usr/include
    /sw/include # Fink
    /opt/local/include # DarwinPorts
    /opt/csw/include # Blastwave
    /opt/include
    /usr/freeware/include
)

IF(CASACORE_INCLUDE_DIR)
    SET(CASACORE_FOUND "YES")
ELSE()
    SET(CASACORE_FOUND "NO")
ENDIF()

#
# Here we call FIND_PACKAGE() on all of the components
#
foreach(_cc_module ${_casacore_modules_to_process})
  string(TOUPPER ${_cc_module} _cc_module_UC)
  FIND_LIBRARY(${_cc_module_UC}_LIBRARIES
    NAMES ${_cc_module}
    PATHS
    $ENV{CASACORE_DIR}/lib
    $ENV{CASACORE_DIR}
  )
  if (NOT ${_cc_module_UC}_LIBRARIES)
    message("Could not find library ${_cc_module}")
    set(CASACORE_FOUND "NO")
  endif()
  MARK_AS_ADVANCED(${_cc_module_UC}_LIBRARIES)
  list(APPEND CASACORE_LIBRARIES ${${_cc_module_UC}_LIBRARIES})
endforeach()

if ("${CASACORE_FOUND}" STREQUAL "NO")
  if (CasaCore_FIND_REQUIRED)
    message(FATAL_ERROR "Casacore not found, please set $CASACORE_DIR")
  else()
    message("Casacore not found, please set $CASACORE_DIR")
  endif()
else()
  message("Found Casacore: ${CASACORE_INCLUDE_DIR}")
endif()

