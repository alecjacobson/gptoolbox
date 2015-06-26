function W = biharmonic_coordinates(V,Ele,b)
  % BIHARMONIC_COORDINATES  Compute the linearly precise generalized
  % coordinates as described in "Linear Subspace Design for Real-Time Shape
  % Deformation" [Wang et al. 2015] **not** to be confused with "Bounded
  % Biharmonic Weights ..." [Jacobson et al. 2011] or "Biharmonic Coordinates"
  % [Weber et al. 12].
  %
  % W = biharmonic_coordinates(V,Ele,b)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   Ele  #Ele by dim+1 list of simplicial element indices into V
  %   b  #b list of indices into V of control points
  % Outputs:
  %   W  #V by #b list of coordinates
  %
  L = cotmatrix(V,Ele);
  [DD,TE] = normal_derivative(V,Ele);
  AT = sparse( ...
    TE(:), ...
    repmat(1:size(TE,1),1,size(TE,2))', ...
    1,size(V,1),size(TE,1));
  [~,C] = on_boundary(Ele);
  DD(~C,:) = 0;
  N = (0.5*AT*DD);
  K = L + N;
  M = massmatrix(V,Ele);
  A = K' * (M\K);
  W = min_quad_with_fixed(A,[],b,eye(numel(b)));
end
