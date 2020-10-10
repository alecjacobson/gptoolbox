function A = quad_facet_adjacency_matrix(Q)
  %
  % A = quad_facet_adjacency_matrix(Q)
  %
  % Inputs:
  %   Q   #Q by 4 list of quad indices
  % Outputs:
  %   A  #Q by #Q adjacency matrix
  %
  allE = [Q(:,1:2);Q(:,2:3);Q(:,3:4);Q(:,[4 1])];
  [E,~,EMAP] = unique(sort(allE,2),'rows');
  E2Q = sparse(EMAP,repmat(1:size(Q,1),1,4)',1,size(E,1),size(Q,1));
  A = (E2Q'*E2Q)>0;
end
