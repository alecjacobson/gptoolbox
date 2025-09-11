> This directory contains mex functions that must be _compiled_ before they can be
called from matlab.

Typically the C++ code for function `myfunc` will be in a file called
`myfunc.cpp` and accompanied by a _blank_ file `myfunc.m`, which simply contains
comments explaining call usage of the function (so that `help myfunc` returns
something useful).

------

**CONTENTS:**
- [Compiling (mexing)](#compiling-mexing)
  - [bash:](#bash)
  - [powershell (or cmd):](#powershell-or-cmd)
  - [Dependencies](#dependencies)
    - [TODO:](#todo)
- [Troubleshooting](#troubleshooting)
  - [MATLAB Versions](#matlab-versions)

------

## Compiling (mexing)

I've abandoned trying to build a "pure-matlab" build system. I now use cmake. So
the build routine is:

### bash:
```bash
    mkdir build 
    cd build
    cmake ..
    make
```

### powershell (or cmd):
```cmd
    mkdir build && cd build
    cmake ..
    cmake --build . --config Release
```
    

This will output the mex functions in this (`mex/`) directory.

CMake's `FindMatlab.cmake` is not very good. You might have to do something like:

    cmake ../ -DMatlab_ROOT_DIR=/apps/matlab-R2019b/

### Dependencies 

Nearly all of the functions depend on stl, Eigen and libigl.  Beyond that some
may depend on CGAL, Embree, and El Topo. The `cmake ..` command above should
take care of _downloading_ these dependencies into `gptoolbox/mex/external/` using vcpkg (including vcpkg itself).

You may already have vcpgk installed in your system. In order to use it, call cmake with additional parameter on configuration step, pointing at your vcpgk executable:
```
    cmake .. -D
```

#### TODO:
+ [ ] As soon as libigl is updated to 2.4.0 in vcpkg (see https://github.com/microsoft/vcpkg/pull/26029), we could switch to vcpkg dependency management also regarding to libigl.

## Troubleshooting

### MATLAB Versions

Cmake unfortunately requires a hardcoded version mapping between MATLAB's numerical
versions (e.g., 9.5) and their named versions (e..g, 2018b). This mapping gets
stale a couple of times a year, so if you find that cmake can't figure out your
MATLAB version, and you have a fairly recent version of MATLAB, it's very likely
that you'll need to update the mapping. To do so, update the `MATLAB_VERSIONS_MAPPING`
variable in `mex/cmake/FindMATLAB.cmake` and add your version number. You can
check your MATLAB version using the `version` command.

If you get an error like

    CMake Error at ... (message):
    Could NOT find Matlab (missing: Matlab_INCLUDE_DIRS Matlab_MEX_LIBRARY
    Matlab_MEX_EXTENSION Matlab_ROOT_DIR MEX_COMPILER MX_LIBRARY ENG_LIBRARY)
    (found version "NOTFOUND")

Then you might consider setting `Matlab_ROOT_DIR` explicitly to help cmake find
MATLAB.

For example, recently I'm using

   rm -f CMakeCache.txt && cmake ../ -DMatlab_ROOT_DIR=/Applications/MATLAB_R2024a.app -DMatlab_MEX_EXTENSION="mexmaca64" -DMATLAB_FIND_DEBUG=ON -DCMAKE_OSX_ARCHITECTURES=arm64

