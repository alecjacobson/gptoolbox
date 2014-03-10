% This is a script file for compiling the mex versions of the Steerable
% Pyramid Tools.
% 
% Usage:>> compilePyrTools
%
% Tested for gcc and lcc.
%
% Rob Young, 9/08

mex upConv.c convolve.c wrap.c edges.c
mex corrDn.c convolve.c wrap.c edges.c
mex histo.c
%mex innerProd.c
mex pointOp.c
mex range2.c
