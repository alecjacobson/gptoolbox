function N = is_vertex_nonmanifold(F)
  % IS_VERTEX_NONMANIFOLD Determine if each vertex in a mesh with faces F is
  % not manifold: incident faces form a single connected component with respect
  % to shared edges.
  %
  % Inputs:
  %   F  #F by 3 list of triangle indices
  % Outputs:
  %   N  #V list of bools whether vertices are not manifold
  %

  function A = face_adjacency(F)
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
    % A(i,j) = 1 means face i and j share an edge
    A = (FE'*FE)>0;
  end

  n = max(F(:));
  m = size(F,1);
  % This is supposedly memory inefficient, but the vectorization in matlab seems
  % to make up for it even in the case of very large meshes m=2,000,000
  vectorized = true;

  if vectorized 
    V2F = sparse(F(:),repmat(1:m,1,3)',1,n,m);
    A = face_adjacency(F);
    [I,J] = find(V2F);
    % UF(i,j) = 1, unique face-for-vertex i is copy of face j
    UF = sparse((1:numel(J))',J,1,numel(J),m);
    
    % U(i,j) = 1, unique face-for-vertex i incident on unique face-for-vertex via
    % single face +1 edge. Problem is that now faces for different vertices may
    % talk to each other.
    U = (UF*A*UF')>0;
    UV = sparse((1:numel(J))',I,1,numel(J),n);
    U = U & (UV*UV')>0;
    % U(i,j) = 1 means unique face-for-vertex i comes from same vertex AND shares
    % edge with face-for-vertex j
    [~,C] = conncomp(U);
    % V2C(i,j) = 1 means vertex i is in component j
    V2C = sparse(I,C,1,n,max(C))>0;
    N = sum(V2C,2)~=1;
  else
    N = zeros(n,1);
    % need transpose for random access
    F2V = sparse(repmat(1:m,1,3)',F(:),1,m,n);
    for v = 1:n
      % face(s) f incident on v
      %f = find(V2F(v,:));
      f = find(F2V(:,v));
      Ff = F(f,:);
      % Need to reindex or memory will explode when building Av
      U = unique(Ff);
      IM = sparse(U,1,1:numel(U));
      Ff = reshape(IM(Ff),size(Ff));
      Av = face_adjacency(Ff);
      N(v) = conncomp(Av)~=1;
    end
  end
end
