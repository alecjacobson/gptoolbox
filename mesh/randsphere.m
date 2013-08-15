function [N] = randsphere(n,varargin)
% RANDSPHERE  randomly sample n points on a sphere
%
% N = randsphere(n)
% N = randsphere(n,'ParameterName',ParameterValue)
%
% Inputs:
%   n  number of points
%   Optional:
%     'Method' followed by method type as string {'trig'}:
%        'naive'  a non-uniform sampling
%        'rejection'  uniform sampling by rejection
%        'trig'  uniform sampling using Archimedes theorem
%        'normal-deviation'  uniform sampling by controling deviation
%        'stratified'  sampling by stratification of lat/lon
%        'lloyd'  uniform(?) sampling by lloyd relaxation of previous method
%          (this method does not scale)
% Outputs:
%   N  n by 3 list of sample positions on sphere
%
% Todo: This function is badly named. Should be something like random_dir
%
  
  method = 'trig';
  
  v = 1;
  while v <= numel(varargin)
    switch varargin{v}
    case 'Method'
      assert((v+1)<=numel(varargin));
      v = v+1;
      method = varargin{v};
    otherwise
      error(sprintf('Unsupported parameter: %s',varargin{v}));
    end
    v = v+1;
  end


  switch method
  case 'naive'
    % uniformly random 3d point, each coord in [-1,1] and normalize
    % Result is non-uniform sampling
    N = normalizerow(2*rand(n,3)-1);
  case 'normal-deviation'
    N = normalizerow(normrnd(zeros(n,3),1));
  case 'rejection'
    N = [];
    % Have to iterate to be sure we get exactly n
    while size(N,1) < n
      % estimate number inliers 
      N_try = rand(ceil((n-size(N,1))/(4*pi/3/8)),3)*2-1;;
      N = [N;normalizerow(N_try(normrow(N_try)<=1,:))];
    end
    N = N(1:n,:);
  case 'trig'
    % "Spherical Sampling by Archimedes’ Theorem" [Shao & Badler 96]
    Z = rand(n,1)*2-1;
    A = rand(n,1)*2*pi;
    R =  sqrt(1-Z.^2);
    N = [Z R.*cos(A) R.*sin(A)];
  case 'stratified'
    % Explicit stratification into floor(sqrt(n)sqrt(n)) ~ n strata
    m = floor(sqrt(n));
    [X,Y] = meshgrid(linspace(0,1-1/(m+1),m));
    % jitter in each strata then uniformly sample any remaining
    J = [[X(:) Y(:)]+1*rand(m^2,2)/(m+1);rand(n-m^2,2)];
    Z = J(:,1)*2-1;
    A = J(:,2)*2*pi;
    R =  sqrt(1-Z.^2);
    N = [Z R.*cos(A) R.*sin(A)];
  case 'lloyd'
    % initialize with sampling
    N = randsphere(n,'Method','stratified');
    lloyd_iter = 10;
    for l = 1:lloyd_iter
      % single lloyd iteration
      F = convhulln(N);
      M = massmatrix(N,F,'voronoi');
      A = adjacency_matrix(F);
      A = A*M;
      %A = bsxfun(@rdivide,A,sum(A,2));
      A = spdiags (1./sum (A,2), 0, size(A,1), size(A,1)) * A ;
      N = A*N;
      % subtract off center of mass  (needed for small n)
      N = bsxfun(@minus,N,diag(M)'*N./sum(diag(M)));
      N = normalizerow(N);
    end 
  otherwise
    error(sprintf('Unsupported Method: %s',method));
  end

end
