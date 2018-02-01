This directory contains mex functions that must be _compiled_ before they can be
called from matlab.

Typically the C++ code for function `myfunc` will be in a file called
`myfunc.cpp` and accompanied by a _blank_ file `myfunc.m`, which simply contains
comments explaining call usage of the function (so that `help myfunc` returns
something useful).

## Compiling (mexing)

I've abandoned trying to build a "pure-matlab" build system. I now use cmake. So
the build routine is (from the bash command line):

    mkdir build
    cd build
    cmake ..
    make 

This will output the mex functions in this (`mex/`) directory.

### Dependencies 

Some of these functions depend on:

 - cmake
 - c++11
   - VS2012 or newer (windows users only)
 - Eigen
 - libigl
 - CGAL
   - boost
 - Embree
 - El Topo
 - Cork
 - Mac OS X Foundation and AppKit frameworks

### libigl

Libigl is by default a _header only_ library. You do not need to compile it to
use it (though you do need to compile and link to any dependencies, e.g. CGAL).

### cmake says Matlab not found?!

cmake has a built in Module to find Matlab. Unfortunately, it only finds
versions of Matlab that have been released at the time you installed cmake. So
if you have a newer Matlab version than your cmake you can:

1. update cmake (e.g., `brew upgrade cmake`), 
2. add your version to the Module: add a line to `/usr/local/share/cmake/Modules/FindMatlab.cmake`, or
3. directly tell cmake where matlab: (e.g., instead of issuing `cmake ..`, issue
   `cmake .. -DMatlab_ROOT_DIR=/Applications/MATLAB_R2018b.app/`)
