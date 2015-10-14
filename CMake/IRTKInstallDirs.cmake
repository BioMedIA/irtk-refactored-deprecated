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

include(GNUInstallDirs)

foreach(dir BINDIR
            LIBDIR
            INCLUDEDIR
            DATADIR
            MANDIR
            DOCDIR)
  set(IRTK_INSTALL_${dir} "${CMAKE_INSTALL_${dir}}" CACHE PATH "")
  set(IRTK_INSTALL_FULL_${dir} "${CMAKE_INSTALL_FULL_${dir}}" CACHE PATH "")
endforeach()

set(IRTK_INSTALL_CMAKEDIR "${IRTK_INSTALL_LIBDIR}/cmake/IRTK" CACHE PATH "")
set(IRTK_INSTALL_FULL_CMAKEDIR "${IRTK_INSTALL_FULL_LIBDIR}/cmake/IRTK" CACHE PATH "")

set(IRTK_INSTALL_EXAMPLEDIR "${IRTK_INSTALL_DOCDIR}/examples" CACHE PATH "")
set(IRTK_INSTALL_FULL_EXAMPLEDIR "${IRTK_INSTALL_FULL_DOCDIR}/examples" CACHE PATH "")
