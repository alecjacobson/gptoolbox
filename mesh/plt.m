function p = plt(V,varargin)
  % PLT
  %
  % p = plt(V,varargin)
  %
  % Shorthand for:
  %  p = plot(V(:,1),V(:,2), ...

  switch size(V,2)
  case 3
    p = plot3(V(:,1),V(:,2),V(:,3),varargin{:});
  case 2
    p = plot(V(:,1),V(:,2),varargin{:});
  end
end 
