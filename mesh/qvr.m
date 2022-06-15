function q = qvr(X,U,varargin)
  % QVR
  %
  % q = qvr(X,U,varargin)
  %
  % Shorthand for:
  %   q = quiver(X(:,1),X(:,2),U(:,1),U(:,2),...)

  if size(X,1) == 1
    X = repmat(X,size(U,1),1);
  end
  if size(U,1) == 1
    U = repmat(U,size(X,1),1);
  end

  switch size(X,2)
  case 3
    qq = quiver3(X(:,1),X(:,2),X(:,3),U(:,1),U(:,2),U(:,3),varargin{:});
  case 2
    qq = quiver(X(:,1),X(:,2),U(:,1),U(:,2),varargin{:});
  end
  
  if nargout>0, q = qq; end
end
