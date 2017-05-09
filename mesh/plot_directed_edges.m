function varargout = plot_directed_edges(V,E,varargin)
  % PLOT_DIRECTED_EDGES Plot direct edges using quiver
  %
  % p = plot_directed_edges(V,E)
  % p = plot_directed_edges(V,E,LINESPEC,...)
  %
  % Inputs:
  %   V  #V by dim list of edge vertiecs
  %   E  #E by 2 list of edges
  % Outputs:
  %   p  plot handle
  %
  % See also: plot_edges
  %
  if isempty(E)
    return;
  end

  switch size(V,2)
  case 2
    q = quiver( ...
      V(E(:,1),1), ...
      V(E(:,1),2), ...
      V(E(:,2),1)-V(E(:,1),1), ...
      V(E(:,2),2)-V(E(:,1),2), ...
      0, ...
      varargin{:});
   otherwise
    error('Unsupported dimension: %d',size(V,2));
   end

  if nargout >= 1
    varargout{1} = p;
  end

end

