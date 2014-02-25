gptoolbox - Geometry Processing Toolbox
=======================================

<https://github.com/alecjacobson/gptoolbox/>

This is a toolbox of useful matlab functions for geometry processing. There are
also tools related to constrainted optimization and image processing. Typically
these are utility functions that are not stand alone applications.

The functions have been organized into folders based on their primary
application:

- dithering/
- geometry/
- imageprocessing/
- images/
- matlabPyrTools/
- matrix/
- mesh/
- meshplot/
- old/
- paths/
- plot/
- utility/

## Installation ##
The vast majority of this code is __straight MATLAB__ (`*.m` files). Thus, only
installing MATLAB and adding the qptoolbox directory and its subdirectories to
your MATLAB path is needed for installation:

    addpath(genpath('/absolute/path/to/gptoolbox'))

To make this change permanent, then issue:

    savepath

There are some mex files, whose documentation for installation are included in
respective `README.md` files. 

To enable tab completion on gptoolbox's IO functions issue:

    add_gptoolbox_tab_completion

this takes a second or two (or 30) and then you'll need to restart MATLAB for
it to take effect.

## License ##
Unless marked otherwise all code is Copyright Alec Jacobson 2014.

We will probably move this code under a MPL2 license in the near future.

## Contact ##
The Geometry Processing Toolbox grew out of Alec Jacobson's private codebase
during his PhD, but has benefited a lot from various collaborators at NYU and
ETH Zurich. Now, the Geometry Processing Toolbox is a group endeavour. If
you're intersted in contributing, please contact Alec Jacobson
(alecjacobson@gmail.com).

## Documentation ##
For now, documentation is limited to a per-function basis. For example, to find
documentation for `cotmatrix` issue:

    help cotmatrix
