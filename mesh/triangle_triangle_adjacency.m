function [Fp, Fi] = triangle_triangle_adjacency(F, varargin)
  % TRIANGLE_TRIANGLE_ADJACENCY Build a face adjacency data structure for a
  % **manifold** triangle mesh. From each face we can find where on its
  % neighboring faces it's incident.
  %
  % [Fp, Fi] = triangle_triangle_adjacency(F)
  % [Fp, Fi] = triangle_triangle_adjacency(F, 'CornerIndexing')
  % [Fp, Fi] = triangle_triangle_adjacency(F, 'CyclicIndexing')
  %
  % This function computes triangle-triangle adjacency for a manifold 
  % triangle mesh. It returns, for each triangle and each edge, the index 
  % of the neighboring triangle and the index of the corresponding edge on 
  % that neighbor.
  %
  % Input:
  %   F  list of face indices, #faces by 3
  % Optional argument (mode):
  %   'CornerIndexing' (default) - Uses corner-based indexing convention: the j-th
  %     neighbor corresponds to the edge opposite the j-th vertex of the triangle.
  %   'CyclicIndexing' - Uses edge-based cyclic indexing convention: the j-th neighbor
  %     corresponds to the j-th edge in the order (1->2, 2->3, 3->1).
  % Output:
  %   Fp  #F by 3 matrix, where Fp(i,j) is the index of the neighboring triangle
  %       adjacent to the j-th edge of triangle i. -1 if the j-th edge is a boundary edge.
  %   Fi  #F by 3 matrix, where Fi(i,j) is the local edge index on the neighboring
  %       triangle corresponding to the j-th edge of triangle i. -1 if the j-th edge is a boundary edge.
  %
  % For example:
  %
  % F = [ 1 3 2;
  %       2 3 4];
  % [Fp, Fi] = triangle_triangle_adjacency(F);
  % Fp =
  %     2    -1    -1
  %    -1    -1     1
  % Fi =
  %     3    -1    -1
  %    -1    -1     1
  % [Fp, Fi] = triangle_triangle_adjacency(F,'CyclicIndexing');
  % Fp =
  %     -1     2    -1
  %      1    -1    -1
  % Fi =
  %     -1     3    -1
  %      1    -1    -1

  % Default mode
  mode = 'CornerIndexing';

  % Parse optional argument
  if ~isempty(varargin)
    if ischar(varargin{1}) || isstring(varargin{1})
      mode = varargin{1};
    else
      error('Optional argument must be a string: ''CornerIndexing'' or ''CyclicIndexing''.');
    end
  end

  if size(F,2) == 3

    % get list of edges (1st edges then 2nd edges then 3rd edges)
    E = [F(:,2) F(:,3); F(:,3) F(:,1); F(:,1) F(:,2)];

    % sparse adjacency matrix where (i,j) = k, k is the index of i-->j in edge
    % list if i-->j AND j-->i exist in edge list. This way is slightly faster
    % than building a proper adjacency list first.
    Ei = sparse(E(:,1),E(:,2),1:size(E,1));
    % adjacency list of edges as if edges represent undirected graph
    unadj = Ei>0;
    % "non-adjacency" list, (i,j) > 0 iff i-->j exists but j-->i does not exist
    nonadj = unadj-(unadj');
    adj = ((unadj + nonadj)==1).*Ei;

    % need to get mapping from sparse ordering to original order
    [~, ~, si] = find(adj);
    % this is slightly faster than sorting the above by ii, which is the same
    [~, ~, v] = find(adj');
    % build map from edges to their corresponding faces
    E2F = repmat(1:size(F,1),1,3)';
    % initialize adjacency map to -1
    Fp = -ones(size(E,1),1);
    Fp(si) = E2F(v);
    Fp = reshape(Fp,size(F));
    % build map from edges to edge positions
    I = reshape(repmat((1:3),size(F,1),1),1,3*size(F,1))';
    % use corresponding edges to find positions of correponding faces
    Fi = -ones(size(E,1),1);
    Fi(si) = I(v);
    Fi = reshape(Fi,size(F));
  else
    error('Not supported yet...');
  end

  % Apply cyclic reordering if requested
  if strcmpi(mode, 'CyclicIndexing')
    Fp = [Fp(:,3), Fp(:,1), Fp(:,2)];
    Fi = [Fi(:,3), Fi(:,1), Fi(:,2)];
    Fi = mod(Fi,3)+1;
  elseif ~strcmpi(mode, 'CornerIndexing')
    error('Unknown indexing mode. Use ''CornerIndexing'' or ''CyclicIndexing''.');
  end
end
