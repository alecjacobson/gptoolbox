function q = qvr(X,U,varargin)
  % QVR
  %
  % q = qvr(X,U,varargin)
  %
  % Shorthand for:
  %   q = quiver(X(:,1),X(:,2),U(:,1),U(:,2),...)

  switch size(X,2)
  case 3
    q = quiver3(X(:,1),X(:,2),X(:,3),U(:,1),U(:,2),U(:,3),varargin{:});
  case 2
    q = quiver(X(:,1),X(:,2),U(:,1),U(:,2),varargin{:});
  end
end
