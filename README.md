libHybridGrains
=======================================
Intro to be finished...

Required Dependencies
---------------------

libHybridGrains requires two dependencies for a minimal build:

* [RapidXML](http://rapidxml.sourceforge.net/): An XML parser to read simulation descriptions.

* [Eigen](http://eigen.tuxfamily.org/): A linear algebra library used internally.

We provide a 'get_dependencies.sh' script to automatically download, verify, and extract the supported versions of these libraries.


Recommended Dependencies
------------------------

We recommend a few dependencies for full featured builds:

* [Qt4](http://qt.digia.com/): A user interface library to provide graphical front ends. Available through most standard package managers.

* [HDF5](https://www.hdfgroup.org/HDF5/): A binary file format for configuration and force output. Available through most standard package managers.

* [Python](https://www.python.org): An interpreted language used for extending libHybridGrains's functionality with plugins. Available standard on most platforms. Note that full libHybridGrains test suite requires the installation of the [numpy](http://www.numpy.org) and [h5py](http://www.h5py.org) Python packages.


Quickstart Guide
----------------

To obtain a minimal demo build that simulates colliding triangle meshes:

1. Install Qt4 and CMake. These packages are available from most standard package managers.

2. Clone this repository and change into the project root:

        git clone https://github.com/peterchencyc/libHybridGrains.git
        cd libHybridGrains

3. From the project root, run the script get_dependencies.sh to download, extract, and verify the required dependencies. Note that this script requires the md5sum utility:

        ./get_dependencies.sh

4. Create a build directory under the project root and change into this directory:

        mkdir build
        cd build

5. Run CMake to create the build system with Qt4 and HDF5 enabled:

        cmake -DUSE_QT4=ON -DUSE_HDF5=ON ..

6. Build libHybridGrains:

        make -j

7. Load the example simulation:

        cd hybridgrains2dnewqt4
        ./hybridgrains2dnewqt4 assets/hourglass/hybrid.xml

8. Click "Simulate" to run the simulation!


Building libHybridGrains
---------------

libHybridGrains uses the CMake build system. libHybridGrains is tested regularly against recent versions of the GCC and LLVM on both Linux and OS X. A minimal libHybridGrains build requires only a C and C++ compiler.

The build system can be configured via the command line by running

    ccmake ..

from the build directory.

Options of note include:

* CMAKE_BUILD_TYPE: General build type that enables various optimization and compiler flags. Options are: Debug, Release, RelWithDebInfo, MinSizeRel, Coverage

* STRICT_BUILD: Enables aggressive warnings and treats warnings as errors. Recommended for development.

* USE_HDF5: Enables state and force output via [HDF5](https://www.hdfgroup.org/HDF5/) files.

* USE_OPENMP: Enables OpenMP support. Highly recommended.

* USE_QT4: Enables support for graphical front ends using [Qt4](http://qt.digia.com/).

* USE_PYTHON: Enables support for embedded Python language scripting. Required for kinematic scripting.

* SANITIZER: Enables support for compiler sanitizer modes. Options are: none, address

Citation
---------------
```
@inproceedings{yue2018hybrid,
  title={Hybrid grains: Adaptive coupling of discrete and continuum simulations of granular media},
  author={Yue, Yonghao and Smith, Breannan and Chen, Peter Yichen and Chantharayukhonthorn, Maytee and Kamrin, Ken and Grinspun, Eitan},
  booktitle={SIGGRAPH Asia 2018 Technical Papers},
  pages={283},
  year={2018},
  organization={ACM}
}
```
