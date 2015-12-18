# Installation instructions

## Dependencies

Required:

 - Any C++11 enabled compiler
 - [CMake](https://www.cmake.org/) version 2.8+
 - [Boost](http://www.boost.org/) version 1.48+

Optional:

 - [Eigen](http://eigen.tuxfamily.org/) version 3
 - [LibFLANN](http://www.cs.ubc.ca/research/flann/)
 - [LibLBFGS](http://www.chokkan.org/software/liblbfgs/)
 - [LibPNG](http://www.libpng.org/pub/png/libpng.html)
 - [NifTICLib](http://niftilib.sourceforge.net/)
 - [VTK](http://www.vtk.org/) version 6
 - MATLAB

## Instructions

### Quick build 

The following snippet performs a quick build and installation of the current release of IRTK using the default options. The [git](https://git-scm.com/) program is required to check the source code out. 

```
git clone https://github.com/BioMedIA/IRTK
cd IRTK/
mkdir build && cd build/
cmake -DCMAKE_INSTALL_PREFIX=../install ..
make && make install
cd ../install/
```

The default installation prefix may be adjusted using the `CMAKE_INSTALL_PREFIX` CMake option. If the latter is not specified, IRTK may be installed under a system-wide path (such as `/usr/local` on Linux), wich may require administrative privileges.

### Advanced build

By default, the common modules and applications of IRTK are built using only the required dependencies. This can be modified by specifying optional build options.

Additional components of IRTK can be built using the `BUILD` family of options. For instance, to build the documentation in addition to the default, use:
```
cmake -DBUILD_DOCUMENTATION=ON
```

Support for optional dependencies can be enabled using the `WITH` family of options. For instance, to build IRTK with VTK support, use:
```
cmake -DWITH_VTK=ON
```

### Installation paths

By default, the installation paths for the different components of IRTK follows the [GNU standard](http://www.gnu.org/prep/standards/html_node/Directory-Variables.html):
 
  - Public headers under `include`.
  - Shared and static libraries under `lib`.
  - Applications under `bin`.
  - Other architecture independent data in their respective subdirectory under `share`.

These paths can be individually overridden using the corresponding `IRTK_INSTALL_*DIR` parameters.

## Known issues

### Build with the system VTK on Ubuntu 14.04

It is not possible to build IRTK from source using the system VTK available on Ubuntu 14.04. The packaged version of VTK (6.0.0) does not work with C++11 projects and the fix for it only landed in subsequent versions. A separate build of VTK should be used instead and the location of its configuration file (`VTKConfig.cmake`) should be specified using the `VTK_DIR` variable:
```
cmake -DWITH_VTK=ON -DVTK_DIR=${VTK_CMAKE_DIRECTORY}
```