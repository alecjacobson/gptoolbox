function l = spline_arc_lengths(P,C,varargin)
  % l = spline_arc_lengths(P,C,...)
  % 
  % Inputs:
  %   P #P by dim list of control points
  %   C #C by 4 list of Bezier curve indices into P
  %   Optional:
  %     'Method' followed by:
  %       'spline_to_poly'  use spline_to_poly to convert to a polygonal
  %       representation and then compute edge lengths 
  %       'cubic_arc_length' use cubic_arc_length to compute arc lengths.
  %       Faster, more accurate. {'cubic_arc_length'}
  %     'Tol' followed by tolerance to pass to spline_to_poly {0.005}
  %     'NumQuadrature' followed by number of quadrature points to pass to
  %       cubic_arc_length. Usually in the range [4,200] {10}
  % Outputs:
  %    
  % See also: spline_to_poly

  method = 'cubic_arc_length';
  tol = 0.005;
  nq = 10;
  if ~isempty(varargin) && isnumeric(varargin{1})
    error('Use new option format');
  end
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Method','Tol','NumQuadrature'}, ...
    {'method','tol','nq'});
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

  switch method
  case 'spline_to_poly'
    [V,E,I] = spline_to_poly(P,C,tol);
    El = edge_lengths(V,E);
    l = accumarray(I,El,[size(C,1) 1]);
  case 'cubic_arc_length'
    l = cubic_arc_length(permute(reshape(P(C',:),4,size(C,1),size(P,2)),[1 3 2]),nq);
  end

end
