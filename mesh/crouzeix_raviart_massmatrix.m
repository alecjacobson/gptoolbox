function [M,E,EMAP] = crouzeix_raviart_massmatrix(V,F)
  % CROUZEIX_RAVIART_MASSMATRIX Compute the Crouzeix-Raviart mass matrix where
  % M(e,e) is just the sum of 1/3 the areas of the triangles on either side of
  % an edge e. For tets, edges are now faccets.
  %
  % See for example "Discrete Quadratic Curvature Energies" [Wardetzky, Bergou,
  % Harmon, Zorin, Grinspun 2007]
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by element-size list of triangle indices
  % Outputs:
  %   M  #E by #E edge-based diagonal mass matrix
  %   E  #E by 2 list of edges
  %
  % See also: edge_laplacian, is_boundary_edge, crouzeix_raviart_cotmatrix,
  %   massmatrix
  %

  switch size(F,2)
  case 3
    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    % Map duplicate edges to first instance
    [E,~,EMAP] = unique(sort(allE,2),'rows');
    TA = doublearea(V,F)/2;
    M = sparse(EMAP,EMAP,repmat(TA/3,3,1),size(E,1), size(E,1));
  case 4
    T = F;
    allF = [ ...
      T(:,2) T(:,4) T(:,3); ...
      T(:,1) T(:,3) T(:,4); ...
      T(:,1) T(:,4) T(:,2); ...
      T(:,1) T(:,2) T(:,3); ...
      ];
    vol = volume(V,T);
    [F,~,FMAP] = unique(sort(allF,2),'rows');
    M = sparse(FMAP,FMAP,repmat(vol/4,4,1),size(F,1), size(F,1));
    E = F;
    EMAP = FMAP;
  end
end
