function [A,uE2F,uE] = facet_adjacency_matrix(F)
  % FACET_ADJACENCY_MATRIX  Adjacency matrix between facets determined by
  % whether two facets share an edge.
  %
  % A = facet_adjacency_matrix(F)
  %
  % Inputs:
  %   F  #F by 3 list of triangles
  % Outputs:
  %   A  #F by #F adjacency matrix 
  %
  % See also: adjacency_matrix

  ss = size(F,2);
  assert(ss==3,'Facets should be triangles');
  % List of all "half"-edges: 3*#F by 2
  allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
  % Sort each row
  sortallE = sort(allE,2);
  % IC(i) tells us where to find sortallE(i,:) in uE: 
  % so that sortallE(i,:) = uE(IC(i),:)
  [uE,~,IC] = unique(sortallE,'rows');
  % uE2F(e,f) = 1 means face f is adjacent to unique edge e
  uE2F = sparse(IC(:),repmat(1:size(F,1),1,ss)',1);
  % kill non-manifold edges
  uE2F(sum(uE2F,2)>2,:) = 0;
  % Face-face Adjacency matrix
  A = uE2F'*uE2F;
  % All ones
  A = A>0;
end
