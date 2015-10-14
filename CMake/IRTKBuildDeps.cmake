# The Image Registration Toolkit (IRTK)
#
# Copyright 2008-2015 Imperial College London
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

find_package(Boost 1.48 REQUIRED COMPONENTS math_c99)
include_directories(${Boost_INCLUDE_DIRS})
link_libraries(${Boost_LIBRARIES})

option(WITH_ZLIB "Build with optional support for ZLIB." OFF)
if(WITH_ZLIB)
  find_package(ZLIB REQUIRED)
  include_directories(${ZLIB_INCLUDE_DIRS})
  link_libraries(${ZLIB_LIBRARIES})
  add_definitions(-DHAS_ZLIB -DHAVE_ZLIB)
endif()

option(WITH_MATLAB "Build with optional support for MATLAB." OFF)
if(WITH_MATLAB)
  find_package(Matlab REQUIRED)
  include_directories(${MATLAB_LIBRARIES})
  add_definitions(-DHAS_MATLAB -DHAVE_MATLAB)
endif()

option(WITH_VTK "Build with optional support for VTK." OFF)
if(WITH_VTK)
  set(VTK_COMPONENTS vtkCommonCore
                     vtkCommonDataModel
                     vtkCommonExecutionModel
                     vtkFiltersCore
                     vtkFiltersFlowPaths
                     vtkFiltersGeneral
                     vtkFiltersGeometry
                     vtkFiltersParallel
                     vtkIOGeometry
                     vtkIOLegacy
                     vtkIOPLY
                     vtkIOXML)
  find_package(VTK COMPONENTS ${VTK_COMPONENTS} NO_MODULE REQUIRED)
  include_directories(${VTK_INCLUDE_DIRS})
  link_libraries(${VTK_LIBRARIES})
  add_definitions(-DHAS_VTK -DHAVE_VTK)
endif()

option(WITH_EIGEN "Build with optional support for Eigen." OFF)
if(WITH_EIGEN)
  find_package(Eigen3 REQUIRED)
  if(DEFINED EIGEN3_USE_FILE)
    include(${EIGEN3_USE_FILE})
  else()
    include_directories(${EIGEN3_INCLUDE_DIR})
  endif()
  add_definitions(-DHAS_EIGEN -DHAVE_EIGEN)
endif()

option(WITH_FLANN "Build with optional support for FLANN." OFF)
if(WITH_FLANN)
  find_package(LibFLANN REQUIRED)
  include_directories(${LIBFLANN_INCLUDE_DIRS})
  link_libraries(${LIBFLANN_LIBRARIES})
  add_definitions(-DHAS_FLANN -DHAVE_FLANN)
endif()

option(WITH_LBFGS "Build with optional support for LBFGS." OFF)
if(WITH_LBFGS)
  find_package(LibLBFGS REQUIRED)
  include_directories(${LIBLBFGS_INCLUDE_DIRS})
  link_libraries(${LIBLBFGS_LIBRARIES})
  add_definitions(-DHAS_LBFGS -DHAVE_LBFGS)
endif()

option(WITH_PNG "Build with optional support for PNG." OFF)
if(WITH_PNG)
  find_package(PNG REQUIRED)
  include_directories(${PNG_INCLUDE_DIRS})
  link_libraries(${PNG_LIBRARIES})
  add_definitions(-DHAS_PNG -DHAVE_PNG)
endif()

option(WITH_NIFTI "Build with optional support for NIFTI." OFF)
if(WITH_NIFTI)
  find_package(NifTICLib REQUIRED)
  include_directories(${NIFTICLIB_INCLUDE_DIRS})
  link_libraries(${NIFTICLIB_LIBRARIES})
  add_definitions(-DHAS_NIFTI -DHAVE_NIFTI)
endif()
