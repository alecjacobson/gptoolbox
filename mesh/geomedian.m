function y = geomedian(X,eta,y)
% y = geomedian(X,eta,y)
%
% Inputs:
%   X  #X by dim list of input points
%   eta  #X list of weights
%   y  initial guess
% Outputs:
%   y  dim-vector output point
%


  if nargin < 2
    eta = [];
  end
  if nargin < 3
    y = median(X);
  end
  X0 = X;
  eta0 = eta;

  % make X unique and accumulate eta if needed
  orig_size = size(X,1);
  [X,~,idx] = unique(X,'rows');
  if size(X,1) < orig_size
    if isempty(eta) 
      eta = accumarray(idx,1);
    else
      eta = accumarray(idx,eta);
    end
  end

  % weiszfeld's algorithm modified in "The multivariate L1-median and associated
  % data depth"
  max_iters = 100;
  for iter = 1:max_iters
    y0 = y;
    V = X-y;
    D = normrow(V);
    keep = D > eps;
    if isempty(eta)
      Q = 1 ./ D;
    else
      Q = eta ./ D;
    end

    y_tilde = sum(Q(keep))^-1 * sum( Q(keep) .* X(keep,:) );
    k = find(~keep,1,'first');
    if isempty(k)
      y = y_tilde;
    else
      if isempty(eta)
        eta_y = 1;
      else
        eta_y = eta(k);
      end
      r_y = normrow(Q(keep,:)'*V(keep,:));
      ratio = eta_y/r_y;
      y = (1-ratio)^-1 * y_tilde + min(1,ratio) * y;
      fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    end
    if norm(y-y0,'inf') < 1e-7
      break;
    end
  end
  if iter == max_iters
    warning('geomedian: max iters reached');
  end
end
