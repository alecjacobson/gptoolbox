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

Nearly all of the functions depend on stl, Eigen and libigl.  Beyond that some
may depend on CGAL, Embree, and El Topo. The `cmake ..` command above should
take care of _downloading_ these dependencies into `gptoolbox/mex/external/`.
