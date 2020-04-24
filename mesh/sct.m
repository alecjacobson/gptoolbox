function s = sct(V,varargin)
  % SCT
  %
  % s = sct(V,varargin)
  %
  % Shorthand for:
  %  s = scatter(V(:,1),V(:,2), ...

  switch size(V,2)
  case 3
    s = scatter3(V(:,1),V(:,2),V(:,3),varargin{:});
  case 2
    s = scatter(V(:,1),V(:,2),varargin{:});
  end
end 
