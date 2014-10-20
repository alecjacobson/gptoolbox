function [A,C] = adjacency_dihedral_angle_matrix(V,F)
  % ADJACENCY_DIHEDRAL_ANGLE_MATRIX Build a matrix A such that A(i,j) = θij the
  % dihedral angle between faces F(i,:) and F(j,:) if they're neighbors (share
  % an edge) or 0 if they're not.
  % 
  % A = adjacency_dihedral_angle_matrix(V,F)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into V. F does *not* need to be
  %     manifold
  % Outputs:
  %   A  #F by #F sparse matrix of signed dihedral angles. All entries are in
  %     [0,2*pi] A-->0 as edge is more convex, A-->pi as edge is more flat,
  %     A-->2*pi as edge is more concave (as view from outside, assuming
  %     counter-clockwise orientation.
  %   C  #F by #F sparse matrix revealing C(i,j) = c that face j is incident on
  %     face i across its corner c
  % 
  % Example:
  %   A = adjacency_dihedral_angle_matrix(V,F);
  %   % unsigned dihedral angles
  %   UA = pi*(A~=0)-abs(pi*(A~=0)-A);
  % 
  %   A = adjacency_dihedral_angle_matrix(V,F);
  %   % Adjacency matrix of nearly coplanar neighbors
  %   UA = pi*(A~=0)-abs(pi*(A~=0)-A);
  %   AF = UA>=(pi-1e-5);
  %   % get connected components
  %   C = components(AF);
  %   % plot with unique colors for each coplanar patch
  %   set(tsurf(F,V),'CData',C);
  %   CM = jet(numel(unique(C)));
  %   colormap(CM(randperm(end),:));
  %

  % all edges "both" directions
  allE = [F(:,[2 3]);F(:,[3 1]);F(:,[2 1])];
  % index for each edge to unique edge
  [E,~,IC] = unique(sort(allE,2),'rows');
  % FE(i,j) = 1 means that face j is incident upon edge i
  % so all FE(:,j) == 1 are neighbors
  FE = sparse( ...
    IC(:), ...
    repmat(1:size(F,1),3,1)', ...
    reshape(repmat(1:3,size(F,1),1),3*size(F,1),1), ...
    size(E,1), ...
    size(F,1));

  % precompute all unit normals
  N = normalizerow(normals(V,F));

  A = sparse(size(F,1),size(F,1));
  C = sparse(size(F,1),size(F,1));
  % Peal off one set of pairs at a time
  % There is a chance this is significantly faster by first transposing FE
  while nnz(FE) > 0
    % Get index of first face per row
    [OM,J] = max(FE,[],2);
    I = 1:size(E,1); % Edge index
    I = I(OM~=0);
    J = J(OM~=0); % First face index
    M = OM(OM~=0);
    % Lookup: L(i) = j  --> edge i reveals face j
    L = sparse(I,1,J,size(E,1),1);
    % remove these from FE
    old = nnz(FE);
    FE = FE - sparse(I,J,M,size(E,1),size(F,1));
    new = nnz(FE);
    assert(new<=old);
    [I2,J2,M2] = find(FE);
    % Pair up with first
    J1 = L(I2);
    % get unit normals
    N1 = N(J1,:);
    N2 = N(J2,:);
    % get unit edge vector
    C1 = mod(OM(I2)+2,3)+1;
    C2 = mod(M2+2,3)+1;
    E1 = sub2ind(size(F),J2,mod(M2+1,3)+1);
    E2 = sub2ind(size(F),J2,mod(M2-1+1,3)+1);
    EV = normalizerow(V(F(E2),:)-V(F(E1),:));
    % Don't need unit normals
    % http://en.wikipedia.org/wiki/Dihedral_angle#Alternative_definitions
    D12 = pi-atan2(dot(cross(N1,N2,2),EV,2),dot(N1,N2,2));
    %D12 = pi-acos(sum(N1.*N2,2));
    % append to A: plus is OK here because if two facets share more than one
    % edge they are coplanar so D12 equals 0 in both cases
    A = A + sparse([J1;J2],[J2;J1],[D12;D12],size(F,1),size(F,1));
    C = C + sparse([J2 J1],[J1 J2],[C2 C1],size(F,1),size(F,1));
  end

end
