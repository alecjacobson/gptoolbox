function [BB,BF] = bounding_box(V,varargin)
  % BOUNDING_BOX  Compute the bounding box of a set of points
  % 
  % [BB,BF] = bounding_box(V)
  % [BB,BF] = bounding_box(V,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by dim list of points
  %   Optional:
  %     'Epsilon' followed by epsilon to offset min and max by {0}
  % Outputs:
  %   BB 2^dim list of boundary box vertices
  %   BF #BF by #dim list of facets
  %

  epsilon = 0;
  params_to_variables = containers.Map( ...
    {'Epsilon'}, ...
    {'epsilon'});
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


  dim = size(V,2);
  minV = min(V)-epsilon;
  maxV = max(V)+epsilon;
  C = combn([0 1],dim) == 1;
  BB = bsxfun(@times,C,minV) + bsxfun(@times,~C,maxV);
  
  % lazy faces
  D = DelaunayTri(BB);
  switch size(D.Triangulation,2)
  case 3
    BF = outline(D.Triangulation);
  case 4
    BF = boundary_faces(D.Triangulation);
  otherwise
    error('Not supported');
  end



end
