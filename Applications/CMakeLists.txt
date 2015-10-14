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

include_directories("../Modules/Common/include")
link_libraries(Common)
include_directories("../Modules/Geometry/include")
link_libraries(Geometry)
include_directories("../Modules/Image/include")
link_libraries(Image)
include_directories("../Modules/Registration/include")
link_libraries(Registration)
include_directories("../Modules/Transformation/include")
link_libraries(Transformation)

if(WITH_VTK)
  include_directories("../Modules/PolyData/include")
  link_libraries(PolyData)
endif()

macro(irtk_add_application _name)
  add_executable(${_name} ${ARGN})
  set_target_properties(${_name} PROPERTIES
                        OUTPUT_NAME irtk-${_name})
  install(TARGETS ${_name}
          DESTINATION ${IRTK_INSTALL_BINDIR})
endmacro()

irtk_add_application(dof2image dof2image.cc)
irtk_add_application(dof2mat dof2mat.cc)
irtk_add_application(dofaverage dofaverage.cc)
irtk_add_application(dofinit dofinit.cc)
irtk_add_application(dofinvert dofinvert.cc)
irtk_add_application(ireg ireg.cc)
irtk_add_application(ptransformation ptransformation.cc)
irtk_add_application(transformation transformation.cc)
