function [facets,E] = merge_planar_patches(V,F,varargin)
  % MERGE_PLANAR_PATCHES
  %
  % [facets] = merge_planar_patches(V,F)
  % [facets,E] = merge_planar_patches(V,F,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices
  %   Optional:
  %     'FacetType' either 'pyramid' or {'tetgen'} whether facets index edges
  %     or vertices
  % Outputs:
  %   facets #planar patches cell of facet outline indices
  %   E  #edges by 2 
  %

  facet_type = 'tetgen';

  v = 1;
  while v <= numel(varargin)
    switch varargin{v}
    case 'FacetType'
      assert((v+1)<=numel(varargin));
      v = v+1;
      facet_type = varargin{v};
    otherwise
      error(['Unsupported parameter: ' varargin{v}]);
    end
    v=v+1;
  end

  switch facet_type
  case 'pyramid'
    error('Not yet supported');
  case 'tetgen'
    A = adjacency_dihedral_angle_matrix(V,F);
    % Adjacency matrix of nearly coplanar neighbors
    AF = A>=(pi-1e-5);
    % get connected components
    C = components(AF);

    % loop over connected components
    facets = {};
    for c = 1:max(C)
      O = outline(F(C==c,:));
      while ~isempty(O)
        facets{end+1} = full(outline_loop(O));
        O = O(~all(ismember(O,facets{end}),2),:);
        holes !!
      end
    end
  end

end
