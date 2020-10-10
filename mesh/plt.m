function p = plt(V,varargin)
  % PLT
  %
  % p = plt(V,varargin)
  %
  % Shorthand for:
  %  p = plot(V(:,1),V(:,2), ...

  switch size(V,2)
  case 3
    pp = plot3(V(:,1),V(:,2),V(:,3),varargin{:});
  case 2
    pp = plot(V(:,1),V(:,2),varargin{:});
  end
  
  if nargout>0, p = pp; end
end 
