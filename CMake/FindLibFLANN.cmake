# Any copyright is dedicated to the Public Domain.
# http://creativecommons.org/publicdomain/zero/1.0/
# Author: Ghislain Antony Vaillant

# - Try to find LibFLANN
# Once done this will define
#  LIBFLANN_FOUND - System has LibFLANN
#  LIBFLANN_INCLUDE_DIRS - The LibFLANN include directories
#  LIBFLANN_LIBRARIES - The libraries needed to use LibFLANN
#  LIBFLANN_DEFNINITIONS - Compiler switches required for using LibFLANN
#
# The content of this file is inspired by the following example:
# https://cmake.org/Wiki/CMake:How_To_Find_Libraries

find_package(PkgConfig)
pkg_check_modules(PC_LIBFLANN QUIET flann)
set(LIBFLANN_DEFINITIONS ${PC_LIBFLANN_CFLAGS_OTHER})

find_path(LIBFLANN_INCLUDE_DIR flann/flann.hpp
          HINTS ${PC_LIBFLANN_INCLUDEDIR}
                ${PC_LIBFLANN_INCLUDE_DIRS})

find_library(LIBFLANN_LIBRARY NAMES flann_cpp
             HINTS ${PC_LIBFLANN_LIBDIR}
                   ${PC_LIBFLANN_LIBRARY_DIRS})

set(LIBFLANN_LIBRARIES ${LIBFLANN_LIBRARY})
set(LIBFLANN_INCLUDE_DIRS ${LIBFLANN_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibFLANN  DEFAULT_MSG
                                  LIBFLANN_LIBRARY LIBFLANN_INCLUDE_DIR)

mark_as_advanced(LIBFLANN_INCLUDE_DIR LIBFLANN_LIBRARY)
