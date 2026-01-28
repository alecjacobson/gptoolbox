function [x0,lp] = chebyshev_center(P)
  % x0 = chebyshev_center(P)
  %
  % Input:
  %   P  #P by dim+1 matrix of halfspaces, where each row is a halfspace
  % Output:
  %   x0  dim-dimensional point that is the Chebyshev center of the convex
  %     polytope

  % norm_vector = np.reshape(np.linalg.norm(halfspaces[:, :-1], axis=1),
  % (halfspaces.shape[0], 1))
  norm_vector = sqrt(sum(P(:,1:end-1).^2,2));
  % c = np.zeros((halfspaces.shape[1],))
  c = zeros(size(P,2),1);
  % c[-1] = -1
  c(end) = -1;
  % A = np.hstack((halfspaces[:, :-1], norm_vector))
  % b = - halfspaces[:, -1:]
  A = [P(:,1:end-1) norm_vector];
  b = -P(:,end);
  % res = linprog(c, A_ub=A, b_ub=b, bounds=(None, None))
  %tic;
  %max_iters = 100;
  %for iter = 1:max_iters
  options = optimoptions( ...
    'linprog','Display','off','Algorithm','dual-simplex-highs' ....
    ... % Improves robustness in zero volume case (from 2^-29 to 2^-30) in a test
    ,'OptimalityTolerance',1e-10,'ConstraintTolerance',1e-10 ...
    );
  [res,fval,exitflag,output] = linprog(c,A,b,[],[],[],[],options);
  %end
  %fprintf('%g secs\n',toc/max_iters);

  %tic;
  %max_iters = 100;
  %for iter = 1:max_iters
  %scipy_linprog(c,A,b);
  %end
  %fprintf('%g secs\n',toc/max_iters);

  x0 = res(1:end-1)';
  %y = res(end);
  %save('prob.mat','A','b','c','res');
  if nargout>=2
    lp = struct('res',res,'fval',fval,'exitflag',exitflag,'output',output);
  end
end

