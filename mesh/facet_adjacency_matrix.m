function [A,uE2F,uE] = facet_adjacency_matrix(F,varargin)
  % FACET_ADJACENCY_MATRIX  Adjacency matrix between facets determined by
  % whether two facets share an edge.
  %
  % A = facet_adjacency_matrix(F)
  % A = facet_adjacency_matrix(F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   F  #F by 3 list of triangles
  %   Optional:
  %     'ManifoldOnly' followed by whether to only consider adjacency across
  %     manifold edges (valence <=2)
  % Outputs:
  %   A  #F by #F adjacency matrix 
  %   uE2F  #E by #F matrix so that (e,f) = 1 means face f is adjacent to
  %     unique edge e
  %   uE  #E by 2 list of unique edges
  %
  % See also: adjacency_matrix

  manifold_only = false;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'ManifoldOnly'},{'manifold_only'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  ss = size(F,2);
  switch ss
  case 3
    % List of all "half"-edges: 3*#F by 2
    allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
    % Sort each row
    sortallE = sort(allE,2);
    % IC(i) tells us where to find sortallE(i,:) in uE: 
    % so that sortallE(i,:) = uE(IC(i),:)
    [uE,~,IC] = unique(sortallE,'rows');
    % uE2F(e,f) = 1 means face f is adjacent to unique edge e
    uE2F = sparse(IC(:),repmat(1:size(F,1),1,ss)',1);
  case 2
    % We're really dealing with edges, so E-->vertices, F-->edges
    uE2F = sparse(F,repmat(1:size(F,1),2,1)',1);
  end

  % kill non-manifold edges
  if manifold_only
    uE2F(sum(uE2F,2)>2,:) = 0;
  end
  % Face-face Adjacency matrix
  A = uE2F'*uE2F;
  % All ones
  A = A>0;
end
