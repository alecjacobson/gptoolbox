function [H] = discrete_mean_curvature(V,F)
  % DISCRETE_MEAN_CURVATURE Compute integrated mean curvature at each vertex.
  %
  % H = discrete_mean_curvature(V,F);
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   H  #V list of integrated mean-curvature values
  %
  % Examples:
  %   H = discrete_mean_curvature(V,F);
  %   M = massmatrix(V,F);
  %   tsurf(F,V,'CData',M\H,'EdgeColor','none',fphong);
  %  
  %
  [A,C] = adjacency_dihedral_angle_matrix(V,F);
  [AI,AJ,AV] = find(A);
  [CI,CJ,CV] = find(C);
  assert(isequal(CI,AI));
  assert(isequal(CJ,AJ));
  l = edge_lengths(V,F);
  % index into F(:)
  opp = sub2ind(size(l),CI,CV);
  inc = [ ...
    sub2ind(size(l),CI,mod(CV+1-1,3)+1) ...
    sub2ind(size(l),CI,mod(CV+2-1,3)+1)];
  lV = l(opp);
  % From Keenan's slide: mean_i = ½ ∑{ij} lij φij
  % Extra 0.5 because each edge is counted twice
  % But why extra-extra 0.5?
  H = full(sparse(F(inc),1,0.5*0.5*0.5*repmat((pi-AV).*l(opp),1,2),size(V,1),1));
end
