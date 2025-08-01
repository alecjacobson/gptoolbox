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
  %   txt(V,[],'BackgroundColor','w'); % num2str((1:size(V,1))') by default

  if numel(varargin)==0 || isempty(varargin{1})
    varargin{1} = num2str((1:size(V,1))');
  end

  if nargin==1
    s = txt(V,num2str((1:size(V,1))'));
    return;
  end

  switch size(V,2)
  case 3
    ss = text(V(:,1),V(:,2),V(:,3),varargin{:});
  case 2
    ss = text(V(:,1),V(:,2),varargin{:});
  end
  
  if nargout>0, s = ss; end
end 

