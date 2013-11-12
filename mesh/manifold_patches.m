function [C,A,uE2F,uE] = manifold_patches(F)
  % Compute connected components of facets connected by manifold edges.
  %
  %[C,A,uE2F,uE] = manifold_patches(F)
  %
  % Known bugs: This will detect a moebius strip as a single patch (manifold,
  % non-orientable) and also non-manfiold, yet orientable patches. 
  %
  % Q: Does this find exactly (manifold || orientable) patches?
  %
  % Inputs:
  %   F  #F by simplex-size list of facets
  % Outputs:
  %   C  #F list of component ids
  %   A  #F by #F adjacency matrix
  %

  % simplex size
  ss = size(F,2);
  assert(ss == 3);

  % List of all "half"-edges: 3*#F by 2
  allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
  % Sort each row
  sortallE = sort(allE,2);
  % IC(i) tells us where to find sortallE(i,:) in uE: 
  % so that sortallE(i,:) = uE(IC(i),:)
  [uE,~,IC] = unique(sortallE,'rows');
  % uE2F(e,f) = 1 means face f is adjacent to unique edge e
  uE2F = sparse(IC(:),repmat(1:size(F,1),1,ss)',1);
  full(uE2F)
  % kill non-manifold edges
  uE2F(sum(uE2F,2)>2,:) = 0;
  % Face-face Adjacency matrix
  A = uE2F'*uE2F;
  % All ones
  A = A>0;
  % Connected components are patches
  %C = components(A); % alternative to graphconncomp from matlab_bgl
  [~,C] = graphconncomp(A);

end
