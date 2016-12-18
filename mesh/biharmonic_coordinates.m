function [W,A,K,M,L,N] = biharmonic_coordinates(V,Ele,b)
  % BIHARMONIC_COORDINATES  Compute the linearly precise generalized
  % coordinates as described in "Linear Subspace Design for Real-Time Shape
  % Deformation" [Wang et al. 2015] **not** to be confused with "Bounded
  % Biharmonic Weights ..." [Jacobson et al. 2011] or "Biharmonic Coordinates"
  % [Weber et al. 12].
  %
  % W = biharmonic_coordinates(V,Ele,b)
  % [W,A,K,M,L,N] = biharmonic_coordinates(V,Ele,b)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   Ele  #Ele by dim+1 list of simplicial element indices into V
  %   b  #b list of indices into V of control points
  % Outputs:
  %   W  #V by #b list of coordinates
  %   A  #V by #V quadratic coefficients matrix
  %   K  #V by #V "Laplacian" so that K = L + N, where L is the "usual"
  %     cotangent Laplacian and N computes normal derivatives at boundary
  %     vertices.
  %

  if size(Ele,2) == 3
    assert(isempty(find_ears(Ele)), ...
      'Mesh should not have ears, preprocess Ele with flip_ears');
  end

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
  if nargin<3 || isempty(b)
    W = [];
  else
    W = min_quad_with_fixed(A,[],b,eye(numel(b)));
  end
end
