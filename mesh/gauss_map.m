function [N,E] = gaussmap(V,F)
  % GAUSSMAP Compute the gaussmap surface of a given mesh (V,F)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertices
  %   F  #F by 3 list of mesh triangle indices
  % Outputs:
  %   N  #F by 3 list of mesh normals (points on gauss map)
  %   E  #E by 2 list of indices into N (spherical edges on gauss map)
  %
  %
  % 
  % Example:
  %   [N,E] = gauss_map(V,F);
  %   subplot(1,2,1);
  %   tsurf(F,V,'FaceVertexCData',0.5*N+0.5)
  %   axis equal;
  %   l = sqrt(sum((N(E(:,1),:) - N(E(:,2),:)).^2,2));
  %   E = E(l>0,:);
  %   [N,I] = remove_unreferenced(N,E);
  %   E = I(E);
  %   [N,E] = spherical_subdivision(N,E,'MaxIter',5,'MinLength',0.1);
  %   subplot(1,2,2);
  %   tsurf([E E(:,1)],N, ...
  %     'FaceVertexCData',0.5*N+0.5,'EdgeColor','interp','LineWidth',2);
  %   axis equal;

  N = normalizerow(normals(V,F));

  % List of all "half"-edges: 3*#F by 2
  allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
  % Sort each row
  sortallE = sort(allE,2);
  % IC(i) tells us where to find sortallE(i,:) in uE:
  % so that sortallE(i,:) = uE(IC(i),:)
  [uE,~,IC] = unique(sortallE,'rows');
  % uE2F(e,f) = 1 means face f is adjacent to unique edge e
  uE2F = sparse(IC(:),repmat(1:size(F,1),1,3)',1);
  % Face-face Adjacency matrix
  A = uE2F'*uE2F;

  [EI,EJ,~] = find(A);
  E = [EI EJ];

end
