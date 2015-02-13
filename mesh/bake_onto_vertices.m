function Z = bake_onto_vertices(S,B,V,F)
  % BAKE_ONTO_VERTICES "Bake" a continuous function S evaulated on a mesh (V,F)
  % at barycentric coordinates B onto vertex values Z. Similar to "Least
  % Squares Vertex Baking" [Kavan et al. 2011]
  %
  % Z = bake_onto_vertices(S,B,V,F)
  %
  % Inputs:
  %   S  #S by 1 list of scalar values
  %   B  #S by #V list of barycentric coordinates such that each row i has
  %     three nonzero entries hoding barycentric corridinates of the point at
  %     which S(i) is evaluated.
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices
  % Outputs:
  %   Z  #V by 1 list of vertex values.
  % 

  w_smooth = 1e-6;
  w_data = 1-w_smooth;

  % number of mesh vertices
  n = size(V,1);

  % Mass matrix
  M = massmatrix(V,F,'voronoi');

  % smoothness term
  L = cotmatrix(V,F);
  % normalize mass so that w_data and w_smooth are invariant
  M = M./sum(M(:));
  Q_smooth = L'*(M\L);
  l_smooth = zeros(n,1);

  % data term
  I = speye(size(S,1))/size(S,1);
  Q_data = B'*(I*B);
  l_data = -2*B'*I*S;

  % combine energy
  Q = w_data*Q_data + w_smooth*Q_smooth;
  l = w_data*l_data + w_smooth*l_smooth;

  %Z = min_quad_with_fixed(Q,l,[],[]);
  lx = zeros(n,1);
  ux = ones(n,1);
  %Z = min_quad_with_fixed_active_set(Q,l,[],[],[],[],[],[],lx,ux);
  Z = quadprog(2*Q,l,[],[],[],[],lx,ux);
  Z(Z>1) = 1;
  Z(Z<0) = 0;
end
