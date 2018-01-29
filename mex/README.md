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
