===========================  matlabPyrTools ============================

This package contains some MatLab tools for multi-scale image
processing.  Briefly, the tools include:
  - Recursive multi-scale image decompositions (pyramids), including
    Laplacian pyramids, QMFs, Wavelets, and steerable pyramids.  These
    operate on 1D or 2D signals of arbitrary dimension.  Data
    structures are compatible with the MatLab wavelet toolbox.
  - Fast 2D convolution routines, with subsampling and boundary-handling.
  - Fast point-operations, histograms, histogram-matching.
  - Fast synthetic image generation: sine gratings, zone plates, fractals, etc.
  - Display routines for images and pyramids.  These include several
    auto-scaling options, rounding to integer zoom factors to avoid 
    resampling artifacts, and useful labeling (dimensions and gray-range).

The package is available as a gnu-zipped UNIX "tar" file, accessible
from the web page:   http://www.cns.nyu.edu/~lcv/software.html

The code was originally written in Matlab version 4.2, and continues
to work in new versions (as of 12/09).  To use the code (these lines
are for UNIX):
  1) gunzip matlabPyrTools.tar.gz  	# unpack g'zipped file
  2) tar tvf matlabPyrTools.tar       	# view contents
  3) tar xvf matlabPyrTools.tar       	# extract into  directory "matlabPyrTools"
  4) rm matlabPyrTools.tar 		# delete tarfile
  5) Run matlab, and execute:
      addpath(<full-pathname-of-matlabPyrTools>);
      help matlabPyrTools

A few functions are actually MEX interfaces to C code.  These are
contained in the subdirectory called MEX.  The MEX files have been
tested on Sun (Solaris), LinuX (on an Intel platform), and Macintosh
OSX (on PowerPC and Intel), but should not be difficult to compile on
most other platforms.  Source code is included in the MEX directory,
as well as Make files.  Pre-compiled versions are included for a
number of platforms.  To compile on your platform, simply run
compilePyrTools.m which is located in the MEX subdirectory.

To make sure these are in your matlab path, you can do *one* of the
following:
  1) Create a symbolic link (or macintosh "alias") for the relavent files 
     in the main matlabPyrTools directory,   or
  2) Copy the relavent files into the main matlabPyrTools directory,  or 
  3) Put the MEX subdirectory in your matlab path: addpath('matlabPyrTools/MEX');

Some example script files showing usage of the code are in the
directory <dir>/TUTORIALS.  There is a README file in that directory
describing the contents.

Incremental changes/updates to the code are documented in the ChangeLog file.

Comments/Suggestions/Bugs to:
  Eero P. Simoncelli
  Center for Neural Science, and
  Courant Institute for Mathematical Sciences
  New York University
  eero.simoncelli@nyu.edu
  http://www.cns.nyu.edu/~eero/
