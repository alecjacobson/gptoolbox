function [C,PI,PC] = dual(V,F)
  % DUAL Construct the dual polygonal mesh of a given triangle mesh
  %
  % [C,PI,PC] = dual(V,F)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of indices into rows of V
  % Outputs:
  %   C  #C by dim list dual vertex positions
  %   PI  #PI stream of polygon indices into rows of C
  %   PC  #V+1 list of cumulative sum of dual face valences
  % 
  % See also: polygons_to_triangles
  %
  assert(~any(on_boundary(F),'all'));
  [~,C] = circumradius(V,F);

  % vertex indices
  I = F(:);
  % only keep unique 
  [I,J] = unique(I);
  % index in faces
  IF = mod(J-1,size(F,1))+1;
  % order in face
  IC = floor((J-1)/size(F,1))+1;

  [Fp,Fi] = triangle_triangle_adjacency(F);

  p1 = @(I) mod(I,3)+1;

  %NF = Fp(sub2ind(size(Fp),IF,p1(IC)));
  %NC = p1(Fi(sub2ind(size(Fi),IF,p1(IC))));
  %clf;
  %hold on;
  %tsurf(F,V,'FaceColor','w',falpha(0.8,1));
  %qvr(V(I,:),C(IF,:)-V(I,:),0,'LineWidth',2);
  %qvr(C(IF,:),C(NF,:)-C(IF,:),0,'LineWidth',2);
  %qvr(C(NF,:),V(F(sub2ind(size(F),NF,NC)),:)-C(NF,:),0,'LineWidth',2);
  %hold off;
  %axis equal;
  %view(3);

  % starting face
  IF0 = IF;
  % ledger of vertex-face pairs in order of observation
  L = [];
  while true
    L = [L;I IF];
    N = sub2ind(size(Fp),IF,p1(IC));
    NF = Fp(N);
    NC = p1(Fi(N));
    %clf;
    %hold on;
    %tsurf(F,V,'FaceColor','w',falpha(0.8,1));
    %qvr(V(I,:),C(IF,:)-V(I,:),0,'LineWidth',2);
    %qvr(C(IF,:),C(NF,:)-C(IF,:),0,'LineWidth',2);
    %qvr(C(NF,:),V(F(sub2ind(size(F),NF,NC)),:)-C(NF,:),0,'LineWidth',2);
    %hold off;
    %axis equal;
    %view(3);
    %pause
    IF = NF;
    IC = NC;

    keep = find(IF ~= IF0);
    IF = IF(keep);
    IC = IC(keep);
    I = I(keep);
    IF0 = IF0(keep);
    if isempty(keep)
      break;
    end
  end

  % stable sort by vertex
  [~,S] = sort(L(:,1));
  L = L(S,:);

  PC = cumsum([0;accumarray(L(:,1),1)]);
  PI = L(:,2);
end
