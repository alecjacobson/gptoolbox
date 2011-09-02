function [O] = outline(F)
  % OUTLINE Find outline (boundary) edges of mesh
  %
  % [O] = outline(F)
  % 
  % Input:
  %  F  #F by 3 face list of indices
  % Output:
  %  O  #O by 2 list of outline edges
  % Example:
  %   [V,F] = create_regular_grid(17,17,0,0);
  %   [O] = outline(F);
  %   % extract unique vertex indices on outline
  %   [u,m,n] = unique(O(:));
  %   % original map O = IM(O)
  %   IM = 1:size(V,1);
  %   IM(O(:)) = n;
  %   % list of vertex positions of outline
  %   OV = V(u,:);
  %   % list of edges in OV 
  %   OE = IM(O);
  %   tsurf(F,V);
  %   hold on;
  %   plot( ...
  %     [OV(OE(:,1),1) OV(OE(:,2),1)]', ...
  %     [OV(OE(:,1),2) OV(OE(:,2),2)]', ...
  %     '-','LineWidth',5);
  %   hold off;
  %

  % Find all edges in mesh, note internal edges are repeated
  E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
  % determine uniqueness of edges
  [u,m,n] = unique(E,'rows');
  % determine counts for each unique edge
  counts = accumarray(n(:), 1);
  % extract edges that only occurred once
  O = u(counts==1,:);
end
