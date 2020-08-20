function s = txt(V,varargin)
  % TXT
  %
  % s = txt(V,varargin)
  %
  % Shorthand for:
  %  s = text(V(:,1),V(:,2), ...
  %
  % Example:
  %   txt(V,num2str((1:size(V,1))'));

  switch size(V,2)
  case 3
    ss = text(V(:,1),V(:,2),V(:,3),varargin{:});
  case 2
    ss = text(V(:,1),V(:,2),varargin{:});
  end
  
  if nargout>0, s = ss; end
end 

