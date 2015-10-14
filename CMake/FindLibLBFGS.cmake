# Any copyright is dedicated to the Public Domain.
# http://creativecommons.org/publicdomain/zero/1.0/
# Author: Ghislain Antony Vaillant

# - Try to find LibLBFGS
# Once done this will define
#  LIBLBFGS_FOUND - System has LibLBFGS
#  LIBLBFGS_INCLUDE_DIRS - The LibLBFGS include directories
#  LIBLBFGS_LIBRARIES - The libraries needed to use LibLBFGS
#  LIBLBFGS_DEFNINITIONS - Compiler switches required for using LibLBFGS
#
# The content of this file is inspired by the following example:
# https://cmake.org/Wiki/CMake:How_To_Find_Libraries

find_package(PkgConfig)
pkg_check_modules(PC_LIBLBFGS QUIET liblbfgs)
set(LIBLBFGS_DEFINITIONS ${PC_LIBLBFGS_CFLAGS_OTHER})

find_path(LIBLBFGS_INCLUDE_DIR lbfgs.h
          HINTS ${PC_LIBLBFGS_INCLUDEDIR}
                ${PC_LIBLBFGS_INCLUDE_DIRS})

find_library(LIBLBFGS_LIBRARY NAMES lbfgs liblbfgs
             HINTS ${PC_LIBLBFGS_LIBDIR}
                   ${PC_LIBLBFGS_LIBRARY_DIRS})

set(LIBLBFGS_LIBRARIES ${LIBLBFGS_LIBRARY})
set(LIBLBFGS_INCLUDE_DIRS ${LIBLBFGS_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(LibLBFGS  DEFAULT_MSG
                                  LIBLBFGS_LIBRARY LIBLBFGS_INCLUDE_DIR)

mark_as_advanced(LIBLBFGS_INCLUDE_DIR LIBLBFGS_LIBRARY)
