function K = cotmatrix3(V,T)
  % COTMATRIX3 computes cotangent matrix for 3D tetmeshes, area/mass terms
  % already cancelled out: laplacian mesh operator 
  %
  % This is distinctly NOT following definition that
  % appears in the appendix of: ``Interactive Topology-aware Surface
  % Reconstruction,'' by Sharf, A. et al
  % http://www.cs.bgu.ac.il/~asharf/Projects/InSuRe/Insure_siggraph_final.pdf
  %
  % Instead it is a purely geometric construction. Find more details in Section
  % 1.1 of "Algorithms and Interfaces for Real-Time Deformation of 2D and 3D
  % shapes" [Jacobson 2013]
  %
  % ND derivation given in "A MONOTONE FINITE ELEMENT SCHEME FOR
  % CONVECTION-DIFFUSION EQUATIONS" [Xu & ZIKATANOV 1999]
  %
  % 3D derivation given in "Aspects of unstructured grids and finite-volume
  % solvers for the Euler and Navier-Stokes equations" [Barth 1992]
  %
  % K = cotmatrix(V,T)
  % Inputs:
  %   V  #V x 3 matrix of vertex coordinates
  %   T  #T x 4  matrix of indices of tetrahedral corners
  % Output:
  %   K  #V x #V matrix of cot weights 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also cotmatrix
  %
  warning('Deprecated. Call cotmatrix directly.');
  K = cotmatrix(V,T);
end
