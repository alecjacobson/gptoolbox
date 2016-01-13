function [C,D] = sample_interior(V,F,n,varargin)
  % SAMPLE_INTERIOR  Sample the interior of a solid bounded by the triangle
  % mesh (V,F)
  %
  % C = sample_interior(V,F,n)
  % C = sample_interior(V,F,n,...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   n  number of samples
  %   Optional:
  %     'MinDist' followed by minimum distance inside mesh, negative implies that
  %       it's OK to be outside by that much.
  % Output:
  %   C  n by 3 list of samples
  % 
  %

  % SAMPLE_BOUNDING_BOX  Smample the interour of the bounding box of V
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   n  number of samples
  % Output:
  %   C  n by 3 list of samples
  %
  function C = sample_bounding_box(V,n)
    C = bsxfun(@plus, bsxfun(@times,rand(n,3),max(V)-min(V)),min(V));
  end

  % default values
  % Map of parameter names to variable names
  min_dist = 0;
  params_to_variables = containers.Map( ...
    {'MinDist'}, ...
    {'min_dist'});
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

  out = 1:n;
  while ~isempty(out)
    C(out,:) = sample_bounding_box(V,numel(out));
    %w_out = winding_number(V,F,C(out,:));
    %out = out(abs(w_out)<0.5);
    % Negative is outside
    D(out) = -signed_distance(C(out,:),V,F);
    out = out(D(out)<min_dist);
  end
end
