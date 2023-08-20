function [c,r] = minimal_bounding_sphere(P)
  % MINIMAL_BOUNDING_SPHERE Compute the minimum enclosing sphere around some
  % points P in d-dimensions. 
  %
  % [c,r] = minimal_bounding_sphere(P)
  %
  % Input:
  %   P  #P by d list of points
  % Outputs:
  %   c  d by 1 center point
  %   r  radius
  %
  % Note: This implementation sets up a conic program with O(n d) variables and
  % O(n) cones.

  % Help reduce the number of auxiliary variables and cones below.
  H = convhull(P);
  P = P(H,:);

  % Skimming the literature it seems like there should be a linear programming
  % formulation, but so far I can only see a conic programming problem:
  % 
  %  min  r
  %  c,r 
  %   s.t. ‖p - c‖ < r
  %
  %  min  r
  %  c,r,d
  %    dᵢ = pᵢ-c
  %    ‖dᵢ‖ ≤ r
  %
  %  min  r
  %  c,r,d
  %    dᵢ+c = pᵢ
  %    ‖dᵢ‖ ≤ r
  %
  %  min  r
  %  c,r,d,k
  %    dᵢ+c = pᵢ
  %    kᵢ = r
  %    ‖dᵢ‖ ≤ kᵢ

  if nargin<2
    R = [];
  else
    R = R(H,:);
  end


  d = size(P,2);
  n = size(P,1);

  [r, res] = mosekopt('symbcon echo(0)');
  prob.c = [repmat(0,d,1);1;zeros(n*d,1);zeros(n,1)]; 
  prob.a = [
    repdiag(sparse(ones(n,1)),d) sparse(n*d,1) speye(n*d,n*d) sparse(n*d,n)
    sparse(n,d) ones(n,1) sparse(n,n*d) -speye(n,n)
    ];
  prob.blc = [P(:);zeros(n,1)];
  prob.buc = [P(:);zeros(n,1)];
  prob.blx = [-inf(d,1);0;-inf(n*d+n,1)];
  prob.bux = inf(d+1+n*d+n,1);
  prob.cones.type = repmat(res.symbcon.MSK_CT_QUAD,n,1);
  prob.cones.sub = [
    d+1+n*d+(1:n)
    d+1+((1:n)+((0:d-1)'*n))
    ];
  prob.cones.sub = prob.cones.sub(:);
  prob.cones.subptr = 1:d+1:n*(d+1);
  [ret,res]=mosekopt('minimize echo(0)',prob);
  c = reshape(res.sol.itr.xx(1:d),1,d);
  r = res.sol.itr.xx(d+1);


end
