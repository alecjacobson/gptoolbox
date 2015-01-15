function W = kharmonic(V,F,b,bc,k)
  % KHARMONIC k-harmonic coordinates, "Harmonic Coordinates for Character
  % Articulation" by Joshi et al, and "An Intuitive Framework for Real-Time
  % Freeform Modeling" section 3 "Precomputed basis functions" by Botsch and
  % Kobbelt
  %
  % W = kharmonic(V,F,b,bc)
  % W = kharmonic(V,F,b,bc,k);
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %  b  list of boundary vertices
  %  bc list of boundary conditions, size(boundary) by # handles matrix of
  %    boundary conditions where bc(:,i) are the boundary conditions for the 
  %    ith handle
  %  Optional:
  %    k  power of laplacian {k=1 %harmonic}
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: boundary_conditions, biharmonic_bounded
  %

  if ~exist('k','var')
    % default to harmonic coordinates
    k = 1;
  end

  % number of vertices
  n = size(V,1);
  % number of handles
  m = size(bc,2);

  % Build discrete laplacian and mass matrices used by all handles' solves
  L = cotmatrix(V,F);
  M = massmatrix(V,F);
  % NORMALIZE MASSMATRIX (THIS IS IMPORTANT!!)
  M = M./max(abs(diag(M)));

  % build k-laplacian, quadratic coefficients
  Q = -L;
  for ii = 2:k
    Q = Q*(M\-L);
  end

  % Minimize W'QW subject to W(b,:) = bc
  tic;
  W = min_quad_with_fixed(Q,zeros(n,1),b,bc);
  toc
end
