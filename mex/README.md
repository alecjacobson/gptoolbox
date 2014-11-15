This directory contains mex functions that must be _compiled_ before they can
be called from matlab.

Typically the C++ code for function `myfunc` will be in a file called
`myfunc.cpp` and accompanied by a _blank_ file `myfunc.m`, which simply
contains comments explaining call usage of the function (so that `help myfunc`
returns something useful).

## Compiling (mexing)

Open `myfunc.cpp` and look for an example `mex` call in the header comments.
Something as simple as:

```matlab
mex -o ray_mesh_intersect ray_mesh_intersect.cpp
```

or a more complicated call like 

```matlab
mex -v -largeArrayDims -DMEX ...
  -I/usr/local/igl/libigl/include ...
  -I/opt/local/include/eigen3 ...
  -I/opt/local/include ...
  -L/opt/local/lib ...
  -lCGAL -lCGAL_Core -lgmp -lmpfr ...
  -lboost_thread-mt -lboost_system-mt ...
  -o signed_distance signed_distance.cpp
```

### Dependencies 

Some of these functions depend on:

 - c++ 11
 - Eigen
 - libigl
 - CGAL
   - boost

### libigl

Libigl is by default a _header only_ library. You do not need to compile it to
use it (though you do need to compile and link to any dependencies, e.g. CGAL).

If you see a mex command _linking_ to libigl libraries (e.g.
`-L/usr/local/igl/libigl -ligl`, etc) you may safely remove these, so long as
you **also remove** any definition of the static library flag:
`-DIGL_STATIC_LIBRARY`.


### OpenMP

OpenMP is _still_ not well supported by clang. I've configured my `mexopts.sh`
to use gcc4.7 (installed via macports) and enabled openmp via the `-fopenmp`
flag. Some of these mex files will utilize this if enabled.

### Windows

Currently, none of these mex functions are developed or even known to be
compilable on Windows. If you succeed, drop me a line. If you fail, I'll try to
help, but I do not have access to a windows machine.
