# Any copyright is dedicated to the Public Domain.
# http://creativecommons.org/publicdomain/zero/1.0/
# Author: Ghislain Antony Vaillant

# - Try to find NifTICLib
# Once done this will define
#  NIFTICLIB_FOUND - System has NifTICLib
#  NIFTICLIB_INCLUDE_DIRS - The NifTICLib include directories
#  NIFTICLIB_LIBRARIES - The libraries needed to use NifTICLib
#
# The content of this file is inspired by the following example:
# https://cmake.org/Wiki/CMake:How_To_Find_Libraries

find_package(PkgConfig)
pkg_check_modules(PC_NIFTICLIB QUIET nifticlib)
set(NIFTICLIB_DEFINITIONS ${PC_NIFTICLIB_CFLAGS_OTHER})

find_path(NIFTICLIB_INCLUDE_DIR nifti1.h
          HINTS ${PC_NIFTICLIB_INCLUDEDIR}
                ${PC_NIFTICLIB_INCLUDE_DIRS}
                "/usr/include/nifti")

find_library(NIFTICLIB_CDF_LIBRARY NAMES nifticdf
             HINTS ${PC_NIFTICLIB_LIBDIR}
                   ${PC_NIFTICLIB_LIBRARY_DIRS})

find_library(NIFTICLIB_IO_LIBRARY NAMES niftiio
             HINTS ${PC_NIFTICLIB_LIBDIR}
                   ${PC_NIFTICLIB_LIBRARY_DIRS})

find_library(NIFTICLIB_ZNZ_LIBRARY NAMES znz
             HINTS ${PC_NIFTICLIB_LIBDIR}
                   ${PC_NIFTICLIB_LIBRARY_DIRS})

set(NIFTICLIB_LIBRARIES ${NIFTICLIB_CDF_LIBRARY}
                        ${NIFTICLIB_IO_LIBRARY}
                        ${NIFTICLIB_ZNZ_LIBRARY})
set(NIFTICLIB_INCLUDE_DIRS ${NIFTICLIB_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NifTICLib DEFAULT_MSG
                                  NIFTICLIB_IO_LIBRARY NIFTICLIB_INCLUDE_DIR)
