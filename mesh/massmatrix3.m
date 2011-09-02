function M = massmatrix3(V,T, type)
  % MASSMATRIX3 mass matrix for the mesh given by V and F
  %
  % M = massmatrix3(V,T, type)
  %
  %
  % Inputs:
  %   V #V x 3 matrix of vertex coordinates
  %   T #T x 4  matrix of indices of tetrahedral corners
  %   type  string containing type of mass matrix to compute
  %     'barycentric': diagonal lumped mass matrix obtained by summing 1/3
  %       of volumes of surrounding tets
  %     Not yet supported:
  %       'full': full mass matrix for p.w. linear fem
  %       'voronoi'
  % Output:
  %   M  #V x #V matrix of cot weights 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: massmatrix
  %

  % vertices must be defined in 3D
  assert(size(V,2)==3);

  % should change code below, so we don't need this transpose
  if(size(T,1) == 4)
    warning('T seems to be 4 by #T, it should be #T by 4');
  end
  T = T';

  if strcmp(type,'full')
    error('full not supported yet...')
  elseif strcmp(type,'barycentric')
    a = V(T(1,:),:);
    b = V(T(2,:),:);
    c = V(T(3,:),:);
    d = V(T(4,:),:);
    % http://en.wikipedia.org/wiki/Tetrahedron#Volume
    % volume for each tetrahedron
    v = repmat(abs(dot((a-d),cross(b-d,c-d,2),2))./6./4,1,4);
    % only diagonal elements
    i = [T(1,:) T(2,:) T(3,:) T(4,:)];
    M = sparse(i,i,v,size(V,1),size(V,1));
  elseif strcmp(type,'voronoi')
    error('voronoi not supported yet...')
  else 
    error('bad mass matrix type')
  end
end
