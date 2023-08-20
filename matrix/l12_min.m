function X = l12_min(A,B,b,bc,C,D)
  % L12_MIN Solve the "L12" minimization problem using second-order cone
  % programming.  L12 is the L2,1 norm applied to the transpose of a matrix
  % (i.e., L12 computes the sum of row-wise vector norms).
  %
  % X = l12_min(A,B,b,bc,C,D)
  %
  % min  ‖ (A X - B)ᵀ ‖_2,1
  %  X
  %
  % min  ‖ A X - B ‖_1,2
  %  X
  %
  % min  ‖ A X - B ‖_1,2
  %  X
  % subject to X(b) = bc
  %
  % min  ‖ A X - B ‖_1,2 + ½ ‖ C X - D ‖_F²
  %  X
  % subject to X(b) = bc
  %
  % Inputs:
  %   A  #A by #X (sparse) matrix
  %   B  #A by #B matrix
  %   b  #b list of indices into rows of X
  %   bc  #b by #B list of fixed values corresponding to b
  %   C  #C by #X (sparse) matrix
  %   D  #C by #B matrix
  % Outputs:
  %   X  #X by #B matrix
  %

  if nargin<3
    b = [];
  end
  if nargin<4
    bc = [];
  end
  if nargin<5
    C = [];
    D = [];
    nq = 0;
  else 
    nq = 1;
  end

  nb = size(B,2);
  na = size(A,1);
  nc = size(C,1);
  n = size(A,2);

  % Make non-empty dimensions consistent
  if isempty(C)
    C = sparse(0,n);
  end
  if isempty(D)
    D = sparse(0,nb);
  end

  prob = struct();
  prob.c = [zeros(n*nb + na*nb,1);ones(na,1);zeros(nc*nb,1);ones(nq,1)];
  prob.a = [
    kroneye(A,nb) -speye(na*nb,na*nb) sparse(na*nb,na + nc*nb + nq); ...
    kroneye(C,nb) sparse(nc*nb,na*nb + na) -speye(nc*nb,nc*nb) sparse(nc*nb,nq); ...
    ];
  prob.blc = [reshape(B',[],1);reshape(D',[],1)];
  prob.buc = prob.blc;
  if ~isempty(b)
    bb = b+(0:nb-1)*n;
    prob.blx = -inf(n*nb+na*nb+na+nc*nb+nq,1);
    prob.bux = inf(size(prob.blx));
    prob.blx(bb(:)) = bc(:);
    prob.bux(bb(:)) = bc(:);
  end
  [~, res] = mosekopt('symbcon echo(0)');
  prob.cones.type = [ ...
    repmat(res.symbcon.MSK_CT_QUAD,1,na) ...
    repmat(res.symbcon.MSK_CT_QUAD,1,nq) ...
    ];
  prob.cones.sub = [ ...
      reshape([n*nb+na*nb+(1:na);reshape(n*nb+(1:na*nb),nb,na)],[],1); ...
      n*nb + na*nb + na + [nc*nb+(1:nq) (1:nc*nb)]' ...
    ];
  prob.cones.subptr = [1:(nb+1):(nb+1)*na na+nb*na+(1:nq)];
  [r,res]=mosekopt('minimize echo(0)',prob);
  X = reshape(res.sol.itr.xx(1:n*nb),n,nb);
end
