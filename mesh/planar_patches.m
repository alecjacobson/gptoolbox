function [C,AF] = planar_patches(V,F,varargin)
  % PLANAR_PATCHES Identify connected patches of co-planar faces
  %
  % [C] = planar_patches(V,F)
  % [C] = planar_patches(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices into V
  %   Optional:
  %     'MinDeltaAngle' followd by minimum change in angle between neighboring
  %     facets to be considered co-planar {pi-1e-5}
  %     'Except' followed by list of faces not to consider during component
  %       analysis
  % Outputs:
  %   C  #F list of patch indices
  %   AF  #F by #F  facet-to-facet adjacency matrix. Faces are consider
  %     connected if their combinatorially connected AND have a dihedral
  %     angle close to pi
  %

  min_delta_angle = pi-1e-5;
  except = [];
  params_to_variables = containers.Map( ...
      {'Except','MinDeltaAngle'}, {'except','min_delta_angle'});
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

  A = adjacency_dihedral_angle_matrix(V,F);
  % Adjacency matrix of nearly coplanar neighbors
  UA = pi*(A~=0)-abs(pi*(A~=0)-A);
  AF = UA>=min_delta_angle;
  AF(except,:) = 0;
  AF(:,except) = 0;
  % get connected components
  [~,C] = conncomp(AF);
  % Force C to be a column vector
  C = reshape(C,numel(C),1);
  %tsurf(F,V,'CData',randcycle(C))
end
