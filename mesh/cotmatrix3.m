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

  if(size(T,1) == 4 && size(T,2) ~=4)
    warning('T seems to be 4 by #T, it should be #T by 4');
  end
  % number of mesh vertices
  n = size(V,1);
  % cotangents of dihedral angles
  C = cotangent(V,T);
  %% TODO: fix cotangent to have better accuracy so this isn't necessary
  %% Zero-out almost zeros to help sparsity
  %C(abs(C)<10*eps) = 0;
  % add to entries
  K = sparse(T(:,[2 3 1 4 4 4]),T(:,[3 1 2 1 2 3]),C,n,n);
  % add in other direction
  K = K + K';
  % diagonal is minus sum of offdiagonal entries
  K = K - diag(sum(K,2));
  %% divide by factor so that regular grid laplacian matches finite-difference
  %% laplacian in interior
  %K = K./(4+2/3*sqrt(3));
  %% multiply by factor so that matches legacy laplacian in sign and
  %% "off-by-factor-of-two-ness"
  %K = K*0.5;
  % flip sign to match cotmatix.m
  if(all(diag(K)>0))
    warning('Flipping sign of cotmatrix3, so that diag is negative');
    K = -K;
  end
end

