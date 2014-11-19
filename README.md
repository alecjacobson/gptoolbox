gptoolbox - Geometry Processing Toolbox
=======================================

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
your MATLAB path is needed for installation:

    addpath(genpath('/absolute/path/to/gptoolbox'))

To make this change permanent, then issue:

    savepath

There are some mex files, whose documentation for installation are included in
respective `mex/README.md` file.

To enable tab completion on gptoolbox's IO functions issue:

    add_gptoolbox_tab_completion

this takes a second or two (or 30) and then you'll need to restart MATLAB for
it to take effect.

### Full installation ###

This strives to be full installation instructions, but will no doubt remain
incomplete for some time. Begin by adding paths as above. 

As stated above, most files are straight matlab and will _just run_ if you have
gptoolbox in your path.

#### Compile `/mex` ####
Most of our mex files will depend on
[libigl](https://github.com/libigl/libigl). The following will assume your
usign a "standard" unix-y install of libigl as a static library.

In MATLAB issue:

    cd mex
    compile_qptoolbox_mex

#### Compile `toolbox_fast_marching` ####

In MATLAB issue:

    cd external/toolbox_fast_marching/
    compile_mex

## Dependencies ##
This depends on MATLAB and its various toolbox extensions. Many functions
should also work with Octave, though this has not been tested.

Functions that rely on `quadprog` have been tested and optimized assuming that
the Mosek toolbox has been installed, but should also work with the `quadprog`
in MATLAB's Optimization Toolbox.

Mex files may have other external dependencies (e.g. CGAL, Eigen, libigl). See
their respective READMEs for more information.

## License ##
Unless marked otherwise, all code is Copyright Alec Jacobson 2014.

We will probably switch to a MPL2 license in the near future.

## Contact ##
The Geometry Processing Toolbox grew out of Alec Jacobson's private codebase
during his PhD, but has benefited a lot from various collaborators at NYU and
ETH Zurich. Now, the Geometry Processing Toolbox is a group endeavour. If
you're intersted in contributing, please contact Alec Jacobson
(alecjacobson@gmail.com) or submit a pull request on github.

## Documentation ##
For now, documentation is limited to a per-function basis. For example, to find
documentation for `cotmatrix` issue:

    help cotmatrix
