# gptoolbox - Geometry Processing Toolbox

![](images/gptoolbox-bounce.gif)

[![Build Status](https://travis-ci.org/alecjacobson/gptoolbox.svg?branch=master)](https://travis-ci.org/alecjacobson/gptoolbox)

<https://github.com/alecjacobson/gptoolbox/>

This is a toolbox of useful matlab functions for geometry processing. There are
also tools related to constrainted optimization and image processing. Typically
these are utility functions that are not stand alone applications.

Here's an incomplete list of cool features this matlab toolbox contains:

- wrappers for TetGen, Triangle, QSlim, meshfix
- mesh smoothing
- mesh clean up (remove duplicates, remove unreferenced)
- geodesic distances on triangle and tetrahedral meshes
- mesh quantities and queries (normals, discrete gaussian curvature, list
  boundary edges, topology, angles, dihedral angles etc.)
- mesh deformation (as-rigid-as-possible (ARAP), moving least-squares, etc.)
- mesh parameterization (harmonic, least squares conformal, ARAP, etc.)
- automatic skinning weight computation (bounded biharmonic weights, bone heat)
- 2D triangle mesh from binary image
- Input/Output for many mesh formats
  (.obj,.off,.stl,.wrl,.ply,.mesh,.node,.ele,.poly,.smf,.bdl,.face)
- discrete differential geometry operators for triangle and tetrahedral meshes
  (cotangent Laplacian, gradient, divergence)
- quadratic programming, active set solver
- scribble-based image colorization, diffusion curves
- exact (un)signed distance field computation for meshes
- constructive solid geometry operations on meshes, booleans
- accelerated point location in triangle and tetrahedral meshes
- image dithering
- deep matlab function dependency

The functions have been organized into folders based on their primary
application:

- external/
- imageprocessing/
- images/
- matrix/
- mesh/
- mex/
- utility/
- wrappers/

## Installation ##
The vast majority of this code is __straight MATLAB__ (`*.m` files). Thus, only
installing MATLAB and adding the qptoolbox directory and its subdirectories to
your MATLAB path is needed for installation. Let's assume you cloned gptoolbox
at `/Users/ajx/Repos/gptoolbox/`, then you could issue:

    addpath(strjoin(strcat(['/Users/ajx/Repos/gptoolbox/'],{'external','imageprocessing', 'images', 'matrix', 'mesh', 'mex', 'quat','utility','wrappers'}),':'))

To make this change permanent, then issue:

    savepath

There are some mex files, whose documentation for installation are included in
respective [mex/README.md](mex/README.md) file.

### Full installation ###

This strives to be full installation instructions, but will no doubt remain
incomplete for some time. Begin by adding paths as above. 

As stated above, most files are straight matlab and will _just run_ if you have
gptoolbox in your path.

#### Compile `/mex` ####

Most of our mex files will depend on [libigl](https://github.com/libigl/libigl).
Some depend on [cgal](https://github.com/CGAL/cgal) and
[embree](https://github.com/embree/embree). We recently switched to a
cmake-based build.

See [mex/README.md](mex/README.md) for details.

#### Compile `external/toolbox_fast_marching/` ####

In MATLAB issue:

    cd external/toolbox_fast_marching/
    compile_mex

## Dependencies ##

gptoolbox depends on MATLAB and some of its toolbox extensions. Many functions
should also work with Octave, though this has not been tested.

Functions that rely on `quadprog` have been tested and optimized assuming that
the Mosek toolbox has been installed, but should also work with the `quadprog`
in MATLAB's Optimization Toolbox.

Mex files may have other external dependencies (e.g. CGAL, Eigen, libigl). See
their respective READMEs for more information. When installing mex libraries,
you may need to modify the files in `wrappers/` (such as `path_to_libigl.m`) so
`gptoolbox` knows where to look.

## Attribution

If you use gptoolbox in your academic projects, please cite the papers we
implement as appropriate. To cite the library in general, you could use this
BibTeX entry:

```bibtex
@misc{gptoolbox,
  title = {{gptoolbox}: Geometry Processing Toolbox},
  author = {Alec Jacobson and others},
  note = {http://github.com/alecjacobson/gptoolbox},
  year = {2024},
}
```

## License ##

This work is dual licensed under MIT and Apache 2.

Please take extra note that the `external/` directory contains code from
elsewhere with potentially more restricted licenses. 

Most of the `mex/` directory includes libigl and Eigen (and some include CGAL),
each with their own licenses.

## Contact ##

The Geometry Processing Toolbox grew out of Alec Jacobson's private codebase
during his PhD, but has benefited a lot from various collaborators at NYU, ETH
Zurich, Columbia, University of Toronto and elsewhere. Now, the Geometry
Processing Toolbox is a group endeavour. If you're intersted in contributing,
please contact Alec Jacobson (alecjacobson@gmail.com) or submit a pull request
on github.

## Documentation ##

For now, documentation is limited to a per-function basis. For example, to find
documentation for `cotmatrix`, in MATLAB issue:

    help cotmatrix
