function [L,E,EMAP] = crouzeix_raviart_cotmatrix(V,F)
  % CROUZEIX_RAVIART_COTMATRIX Compute the Crouzeix-Raviart cotangent
  % Laplacian matrix where we use test functions define at edge midpoints.
  %
  % See for example "Discrete Quadratic Curvature Energies" [Wardetzky, Bergou,
  % Harmon, Zorin, Grinspun 2007]
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by element-size list of triangle indices
  % Outputs:
  %   L  #E by #E edge-based sparse cotangent matrix
  %   E  #E by 2 list of edges
  %
  % Examples:
  %   % mesh in (V,F)
  %   [~,~,Ucr] = svd(full(Lcr));
  %   % Convert between edge values and vertex values
  %   E2V = sparse(E(:),repmat(1:size(E,1),1,2)',1,size(V,1),size(E,1));
  %   E2V = bsxfun(@rdivide,E2V,sum(E2V,2));
  %   V2E = sparse(E(:),repmat(1:size(E,1),1,2)',1,size(V,1),size(E,1))';
  %   V2E = bsxfun(@rdivide,V2E,sum(V2E,2));
  %   tsurf(F,[V(:,1:2) E2V*Ucr(:,end-1)])
  %
  %
  % See also: edge_laplacian, is_boundary_edge, crouzeix_raviart_massmatrix,
  %   cotmatrix
  %

  %% check for non-manifold edges
  %S = statistics(V,F,'Fast',true);
  %if S.num_nonmanifold_edges > 0
  %  error(sprintf('There are %d non-manifold edges',num_nonmanifold_edges));
  %end

  % number of vertices
  n = size(V,1);
  % number of elements
  m = size(F,1);
  % simplex size
  ss = size(F,2);
  switch ss
  case 3
    % Compute cotangents: seems to be 0.5*C
    C = 2*cotangent(V,F);

    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    % Map each face-edge to a unique edge
    F2E = reshape(1:3*m,m,3);
    % Lij = -2 cot aij
    %
    % o
    % |\
    % | \
    % |  \
    % |   \
    % i    \
    % |     \
    % |      \
    % |αij    \
    % o----j---o
    %
    %
    LI =  F2E(:, [1 2 3 2 3 1        1 2 3 2 3 1]);
    LJ =  F2E(:, [2 3 1 1 2 3        1 2 3 2 3 1]);
    LV = 2*[-C(:,[3 1 2 3 1 2]) C(:,[3 1 2 3 1 2])];

    % Map duplicate edges to first instance
    [E,~,EMAP] = unique(sort(allE,2),'rows');

    assert(all(size(LI)==size(LJ)));
    assert(all(size(LI)==size(LV)));
    L = sparse(EMAP(LI),EMAP(LJ),LV,size(E,1),size(E,1));
  case 4
    % tets
    T = F;
    C = -2*cotangent(V,T);
    allF = [ ...
      T(:,2) T(:,4) T(:,3); ...
      T(:,1) T(:,3) T(:,4); ...
      T(:,1) T(:,4) T(:,2); ...
      T(:,1) T(:,2) T(:,3); ...
      ];
    % Map each element-face to a unique face
    T2F = reshape(1:4*m,m,4);
    % Lij = -2 lij cot αij
    LI = T2F(:,[1 4 4 4 2 3 2 1 2 3 3 1        1 4 4 4 2 3 2 1 2 3 3 1]);
    LJ = T2F(:,[2 1 2 3 3 1 1 4 4 4 2 3        1 4 4 4 2 3 2 1 2 3 3 1]); 
    LV = [-C(:,[3 4 5 6 1 2 3 4 5 6 1 2]) C(:,[3 4 5 6 1 2 3 4 5 6 1 2])];
    % Map duplicate facets to first instance
    [F,~,FMAP] = unique(sort(allF,2),'rows');

    assert(all(size(LI)==size(LJ)));
    assert(all(size(LI)==size(LV)));
    L = sparse(FMAP(LI),FMAP(LJ),LV,size(F,1),size(F,1));
    E = F;
    EMAP = FMAP;

  otherwise
    error(['Simplex size ' num2str(ss) ' unsupported']);
  end

end
