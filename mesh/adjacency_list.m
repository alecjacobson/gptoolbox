function A = adjacency_list(F,b)
  % ADJACENCY_LIST build adjacency lists for vertices of triangle mesh
  %
  % A = adjacency_list(F)
  % A = adjacency_list(F,b)
  % 
  % Inputs:
  %  F  #F by 3 list of triangle indices
  %  Optional:
  %    b  #b list of indices of vertices for which to build adjancency lists,
  %      default is 1:max(F(:))
  % Output:
  %  A  #b cell array of adjacency lists so that A{i} are the neighbors of b(i)
  %


  indices = 1:max(F(:));

  % if b is not given build for all
  if ~exist('b','var')
    b = indices;
  end

  assert(numel(b) == prod(size(b)));
  % number of vertices for which to build adjacency lists
  m = numel(b);

  % build adjacency matrix
  M = adjacency_matrix(F);

  A = cell(m,1);
  for ii = 1:m
    A{ii} = indices(logical(M(b(ii),:)));
  end


end
