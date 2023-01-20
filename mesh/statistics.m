function S = statistics(V,F,varargin)
  % STATISTICS  Determine a number of statistics about a mesh
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertices
  %   F  #F by 3 list of triangle indices
  %   Optional:
  %     'Fast'  followed by bool. Whether to skip certain expensive statistics
  %     'MinArea'  Followed by factor of bounding_box_diagonal^2 area value for
  %       determining small triangles {1e-7}*bounding_box_diagonal^2
  %     'MinDist'  Followed by factor of bounding box diagonal dist value for
  %       determining close vertices {1e-7}*bounding_box_diagonal
  %     'MinAngle'  Followed by dist value for determining small angles
  %       {0.01 radians}
  % Outputs:
  %   S  Struct containing:
  %     num_faces  Number of triangles in F
  %     num_vertices  Number of vertices in V
  %     num_edges  Number of unique undirected edges
  %     num_small_triangles  Number of triangles with area < min_area
  %     num_small_angles  Number of corners with internal angles < min_angle
  %     num_close_vertices  Number of vertices minus number of unique vertices
  %       upto rounding coordinates by min_dist (~L1 distance)
  %     num_connected_components Number of connected compenents of mesh as
  %       graph (includes each unreferenced vertex as a singlton component)
  %     num_handles  Genus of surface, (2-X)/2, where X is Euler characteristic
  %     euler_characteristic  Euler characteristic
  %     num_boundary_loops  Number of connected components in undirected graph
  %       of boundary edges
  %     num_boundary_edges  Number of boundary edges
  %     num_ears  Number of topological ears (triangles with single internal
  %       edge)
  %     num_nonmanifold_edges Number of non-manifold edges (edges with valence
  %       >2)
  %     num_conflictingly_oriented_edges Number of edges with and even number
  %       of incident faces having conflicting orientation
  %     num_duplicate_vertices  Number vertices minus number of geometrically
  %       unique vertices
  %     num_nonmanifold_vertices Number of non-manifold vertices (vertices
  %       whose incident faces do not form exaclty one connected component)
  %     num_unreferenced_vertices Number of vertices not appearing in F
  %     num_combinatorially_duplicate_faces  Number of faces minus number of
  %       combinatorially unique faces
  %     num_geometrically_degenerate_faces  Number of faces with area == 0.
  %       Note that combinatorially degenerate faces will be included.
  %     num_combinatorially_degenerate_faces  Number of faces with 2 or more of
  %       the same vertices as corners
  %     Not 'Fast'
  %       num_selfintersecting_pairs  Number of self intersecting pairs of
  %         triangles
  %       num_intracomponent_selfintersecting_pairs Number of self intersecting
  %          pairs of triangles _in the same connected component_
  %
  % Examples:
  %   S = statistics(V,F,'Fast',true)
  %   subplot(3,1,1);
  %   hist(doublearea(V,F)/2,100);
  %   title('areas');
  %   subplot(3,1,2);
  %   hist(reshape(internalangles(V,F),[],1),100);
  %   title('internal angles');
  %   subplot(3,1,3);
  %   bar( ...
  %    sort(sparse( ...
  %      connected_components(F),1,diag(massmatrix(V,F,'barycentric'))), ...
  %      'descend'));
  %   title('CC areas');

  fast = false;
  min_area = 1e-7;
  min_dist = 1e-7;
  min_angle = 0.01;

  ii = 1;
  while ii<=numel(varargin)
    switch varargin{ii}
    case 'MinArea'
      assert(ii+1<=numel(varargin));
      ii = ii+1;
      min_area = varargin{ii};
    case 'MinDist'
      assert(ii+1<=numel(varargin));
      ii = ii+1;
      min_dist= varargin{ii};
    case 'MinAngle'
      assert(ii+1<=numel(varargin));
      ii = ii+1;
      min_angle = varargin{ii};
    case 'Fast'
      assert(ii+1<=numel(varargin));
      ii = ii+1;
      fast = varargin{ii};
    otherwise
      error(['Unsupported parameter `' varargin{ii} '`']);
    end
    ii = ii + 1;
  end

  bbd = sqrt(sum((max(V)-min(V)).^2,2));
  
  % To achieve this ordering (or almost this order?):
  % Basic counts
  S.num_faces = 0;
  S.num_vertices = 0;
  S.num_edges = 0;
  % Small things
  S.num_small_triangles = 0;
  S.num_small_angles = 0;
  S.num_kahan_length_violations = 0;
  S.num_close_vertices = 0;
  % Topology related
  S.num_connected_components = 0;
  S.num_handles = 0;
  S.euler_characteristic = 0;
  S.num_boundary_loops = 0;
  S.num_boundary_edges = 0;
  S.num_ears = 0;
  S.num_nonmanifold_edges = 0;
  S.num_conflictingly_oriented_edges = 0;
  S.num_duplicate_vertices = 0;
  S.num_nonmanifold_vertices = 0;
  S.num_unreferenced_vertices = 0;
  % Degeneracy related
  S.num_combinatorially_duplicate_faces = 0;
  S.num_geometrically_degenerate_faces = 0;
  S.num_combinatorially_degenerate_faces = 0;
  if ~fast
    S.num_intracomponent_selfintersecting_pairs = 0;
    S.num_selfintersecting_pairs = 0;
  end

  % easy
  S.num_faces = size(F,1);
  S.num_vertices = size(V,1);
  S.num_edges = size(edges(F),1);
  S.num_combinatorially_duplicate_faces = ...
    S.num_faces - size(unique(sort(F,2),'rows'),1);
  S.num_duplicate_vertices = ...
    S.num_vertices - size(remove_duplicate_vertices(V,0),1);

  C = zeros(size(V,1),1);
  C(1:max(F(:)),:) = connected_components(F);
  C(C==0) = max(C)+(1:sum(C==0));
  is_unreferenced = sparse(F,1,1,size(V,1),1)==0;
  %S.num_unreferenced_vertices = sum(sparse(C,1,1)==1)
  S.num_unreferenced_vertices = sum(is_unreferenced);
  S.num_connected_components = max(C);

  E = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
  % Direct all edges so sortE(:,1) < sortE(:,2)
  sortE = sort(E,2);
  % Adjacency matrix for these "redirected" edges
  DA = sparse(sortE(:,1),sortE(:,2),1,size(V,1),size(V,1));
  % If edge occurs more than once it's non-manifold
  S.num_nonmanifold_edges = nnz(DA>2);
  % If edge only occurs once then it's a boundary
  S.num_boundary_edges = nnz(DA==1);
  if S.num_nonmanifold_edges == 0
    S.num_ears = size(find_ears(F),1);
  else
    rmfield(S,'num_ears');
  end
  % Same adjacency matrix but count -1 if E(:,1)<E(:,2) for and +1 otherwise
  OA = sparse(sortE(:,1),sortE(:,2),1-2*(E(:,1)<E(:,2)),size(V,1),size(V,1));
  % Don't count boundary edges (where "redirected" edge only occured once).
  OA(DA==1) = 0;
  % non-zero count means conflictingly oriented. Note, can get +2-1=0, but this
  % could feasibly be oriented (a meshed self-intersection).
  S.num_conflictingly_oriented_edges = nnz(OA);

  isnmv = is_vertex_nonmanifold(F);
  if numel(isnmv)<size(V,1);
    isnmv(numel(isnmv)+1:size(V,1)) = 0;
  end
  S.num_nonmanifold_vertices = sum(isnmv & ~is_unreferenced);

  [~,BC] = conncomp( (DA==1)+(DA==1)' );
  S.num_boundary_loops = sum(sparse(BC,1,1)>1);

  S.euler_characteristic = S.num_vertices - S.num_edges + S.num_faces;
  b = S.num_boundary_loops;
  % http://en.wikipedia.org/wiki/Genus_(mathematics)#Orientable_surface
  % X = 2 - 2g - b
  % X + b - 2 = - 2g
  % 2-b-X  = 2g
  % g = (2-b-X)/2
  % http://sketchesoftopology.wordpress.com/2008/02/04/genus-euler-characteristic-boundary-components/
  S.num_handles = ...
    (2* ...
      (S.num_connected_components-S.num_unreferenced_vertices) - ...
    S.num_boundary_loops - ...
    (S.euler_characteristic-S.num_unreferenced_vertices))/2 ;

  dblA = doublearea(V,F);
  S.num_geometrically_degenerate_faces = sum(dblA==0);
  S.num_combinatorially_degenerate_faces = ...
    sum((F(:,1)==F(:,2)) | (F(:,2)==F(:,3)) | (F(:,3)==F(:,1)));

  S.num_small_triangles = sum(dblA<2*min_area*bbd*bbd);
  S.num_small_angles = sum(sum(internalangles(V,F)<min_angle));
  sl = sort(edge_lengths(V,F),2,'descend');
  S.num_kahan_length_violations = sum((sl(:,3)-(sl(:,1)-sl(:,2)))<0);

  S.num_close_vertices = ...
    S.num_vertices - size(remove_duplicate_vertices(V,min_dist*bbd),1);

  if ~fast
    V3 = V;
    V3(:,end+1:3) = 0;
    nd = dblA>0;
    Fnd = F(nd,:);
    [~,~,IF] = selfintersect(V3,Fnd,'DetectOnly',true);
    S.num_selfintersecting_pairs = size(IF,1);
    CF = C(F(:,1));
    CFnd = CF(nd)';
    S.num_intracomponent_selfintersecting_pairs =  ...
      sum(CFnd(IF(:,1))== CFnd(IF(:,2)));
  end


end
