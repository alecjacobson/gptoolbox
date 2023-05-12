function [V,F,I,C] = rhombicosidodecahedron()
  %
  % Example:
  %   [V,F,I,C] = rhombicosidodecahedron();
  %   tsurf(F,V,'CData',I,'EdgeColor','none',falpha(1,1));
  %   hold on;
  %   E = sharp_edges(V,F);
  %   tsurf(E,V,'LineWidth',2,'FaceColor','none','EdgeColor','k');
  %   qvr(zeros(size(C)),C*3);
  %   hold off;
  %   axis equal;

  V = [-1/2, -1/2, -1 - sqrt(5)/2;-1/2, -1/2, (2 + sqrt(5))/2;-1/2, 1/2, -1 - sqrt(5)/2;-1/2, 1/2, (2 + sqrt(5))/2;-1/2, -1 - sqrt(5)/2, -1/2;-1/2, -1 - sqrt(5)/2, 1/2;-1/2, (2 + sqrt(5))/2, -1/2;-1/2, (2 + sqrt(5))/2, 1/2;0, (-3 + sqrt(5))^(-1), (-5 - sqrt(5))/4;0, (-3 + sqrt(5))^(-1), (5 + sqrt(5))/4;0, (3 + sqrt(5))/4, (-5 - sqrt(5))/4;0, (3 + sqrt(5))/4, (5 + sqrt(5))/4;1/2, -1/2, -1 - sqrt(5)/2;1/2, -1/2, (2 + sqrt(5))/2;1/2, 1/2, -1 - sqrt(5)/2;1/2, 1/2, (2 + sqrt(5))/2;1/2, -1 - sqrt(5)/2, -1/2;1/2, -1 - sqrt(5)/2, 1/2;1/2, (2 + sqrt(5))/2, -1/2;1/2, (2 + sqrt(5))/2, 1/2;(-5 - sqrt(5))/4, 0, (-3 + sqrt(5))^(-1);(-5 - sqrt(5))/4, 0, (3 + sqrt(5))/4;(-1 - sqrt(5))/4, (-1 - sqrt(5))/2, (-3 + sqrt(5))^(-1);(-1 - sqrt(5))/4, (-1 - sqrt(5))/2, (3 + sqrt(5))/4;(-1 - sqrt(5))/4, (1 + sqrt(5))/2, (-3 + sqrt(5))^(-1);(-1 - sqrt(5))/4, (1 + sqrt(5))/2, (3 + sqrt(5))/4;(-1 - sqrt(5))/2, (-3 + sqrt(5))^(-1), (-1 - sqrt(5))/4;(-1 - sqrt(5))/2, (-3 + sqrt(5))^(-1), (1 + sqrt(5))/4;(-1 - sqrt(5))/2, (3 + sqrt(5))/4, (-1 - sqrt(5))/4;(-1 - sqrt(5))/2, (3 + sqrt(5))/4, (1 + sqrt(5))/4;-1 - sqrt(5)/2, -1/2, -1/2;-1 - sqrt(5)/2, -1/2, 1/2;-1 - sqrt(5)/2, 1/2, -1/2;-1 - sqrt(5)/2, 1/2, 1/2;(-3 + sqrt(5))^(-1), (-5 - sqrt(5))/4, 0;(-3 + sqrt(5))^(-1), (-1 - sqrt(5))/4, (-1 - sqrt(5))/2;(-3 + sqrt(5))^(-1), (-1 - sqrt(5))/4, (1 + sqrt(5))/2;(-3 + sqrt(5))^(-1), (1 + sqrt(5))/4, (-1 - sqrt(5))/2;(-3 + sqrt(5))^(-1), (1 + sqrt(5))/4, (1 + sqrt(5))/2;(-3 + sqrt(5))^(-1), (5 + sqrt(5))/4, 0;(1 + sqrt(5))/4, (-1 - sqrt(5))/2, (-3 + sqrt(5))^(-1);(1 + sqrt(5))/4, (-1 - sqrt(5))/2, (3 + sqrt(5))/4;(1 + sqrt(5))/4, (1 + sqrt(5))/2, (-3 + sqrt(5))^(-1);(1 + sqrt(5))/4, (1 + sqrt(5))/2, (3 + sqrt(5))/4;(1 + sqrt(5))/2, (-3 + sqrt(5))^(-1), (-1 - sqrt(5))/4;(1 + sqrt(5))/2, (-3 + sqrt(5))^(-1), (1 + sqrt(5))/4;(1 + sqrt(5))/2, (3 + sqrt(5))/4, (-1 - sqrt(5))/4;(1 + sqrt(5))/2, (3 + sqrt(5))/4, (1 + sqrt(5))/4;(2 + sqrt(5))/2, -1/2, -1/2;(2 + sqrt(5))/2, -1/2, 1/2;(2 + sqrt(5))/2, 1/2, -1/2;(2 + sqrt(5))/2, 1/2, 1/2;(3 + sqrt(5))/4, (-5 - sqrt(5))/4, 0;(3 + sqrt(5))/4, (-1 - sqrt(5))/4, (-1 - sqrt(5))/2;(3 + sqrt(5))/4, (-1 - sqrt(5))/4, (1 + sqrt(5))/2;(3 + sqrt(5))/4, (1 + sqrt(5))/4, (-1 - sqrt(5))/2;(3 + sqrt(5))/4, (1 + sqrt(5))/4, (1 + sqrt(5))/2;(3 + sqrt(5))/4, (5 + sqrt(5))/4, 0;(5 + sqrt(5))/4, 0, (-3 + sqrt(5))^(-1);(5 + sqrt(5))/4, 0, (3 + sqrt(5))/4];
  F = convhull(V);

  N = normalizerow(normals(V,F));
  %[I,C] = kmeans(N,62);
  [U,~,I] = unique(round(N*2),'rows');
  C = full_sparse([I I I],repmat(1:3,size(I,1),1),N,max(I),3)./accumarray(I,1);
end
