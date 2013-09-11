# Locate casa-core
# This module defines
# WCS_LIBRARIES
# WCS_FOUND, if false, do not try to link to libwcs 
# WCS_INCLUDE_DIR, where to find the headers
#
# $WCS_DIR is an environment variable that would
# correspond to the ./configure --prefix=$WCS_DIR
# used in building libwcs.

FIND_PATH(WCS_INCLUDE_DIR wcslib/wcs.h
    NO_DEFAULT_PATH
    PATHS
    $ENV{WCS_ROOT_DIR}/include
    $ENV{WCS_ROOT_DIR}
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

FIND_LIBRARY(WCS_LIBRARIES
  NAMES wcs
  PATHS
  $ENV{WCS_ROOT_DIR}/lib
  $ENV{WCS_ROOT_DIR}
  )

message("Include directory: ${WCS_INCLUDE_DIR}")
message("Library directory: ${WCS_LIBRARIES}")
if (WCS_INCLUDE_DIR AND WCS_LIBRARIES)
  message("Found Wcs: ${WCS_INCLUDE_DIR}")
else()
  if (Wcs_FIND_REQUIRED)
    message(FATAL_ERROR "Wcs not found, please set $WCS_ROOT_DIR")
  else()
    message("Wcs not found, please set $WCS_ROOT_DIR")
  endif()
endif()

