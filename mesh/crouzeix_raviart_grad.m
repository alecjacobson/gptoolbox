function [G,uE,EMAP] = crouzeix_raviart_grad(V,F)
  % CROUZEIX_RAVIART_GRAD Computes the Gradient matrix for non-conforming
  % Crouzeix-Raviart elements defined over a triangle mesh.
  %
  % Inputs:
  %   V  #V by dim list of triangle vertex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   G  #F*dim by #uE matrix computing piecewise constant gradient vectors at
  %     faces given piecewise linear functions defined by values at edge
  %     mid-points
  %   uE  #uE by 2 list of unique (undirected) edges
  %   EMAP  #F*3 list of indices into uE mapping directed half-edges to unique
  %     edges.
  %
  % Example:
  %   [G,E,EMAP] = crouzeix_raviart_grad(V,F);
  %   L = crouzeix_raviart_cotmatrix(V,F);
  %   GMG = G'*repdiag(diag(sparse(dblA)/2),3)*G;
  %   assert(max(max(abs(L-GMG))) < 1e-12)
  %

  assert(size(F,2) == 3,'Only triangles are supported');
  
  dim = size(V,2);
  V(:,end+1:3) = 0;
  m = size(F,1);
  E = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
  % Map each face-edge to a unique edge
  F2E = reshape(1:3*m,m,3);
  % Map duplicate edges to first instance
  [uE,~,EMAP] = unique(sort(E,2),'rows');

  % Edge vectors
  EV = V(E(:,2),:)-V(E(:,1),:);

  N = normals(V,F);
  dblA = normrow(N);
  U = normalizerow(N);

  % rotate each vector 90 degrees around normal and divide by area
  EVp = -cross(repmat(bsxfun(@rdivide,U,dblA/2),3,1),EV);
  EVp = EVp(:,1:dim);

  G = sparse( ...
    bsxfun(@plus,repmat(1:m,1,3)',m*(0:dim-1)), ...
    repmat(EMAP,1,dim), ...
    EVp+eps, ...
    dim*m,size(uE,1));
    

end
