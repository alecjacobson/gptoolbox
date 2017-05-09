function [V,F] = cookie_cutter(P,E)
  % COOKIE_CUTTER Generate a cookie cutter model from a curve
  %
  % [V,F,W,G] = cookie_cutter(P,E)
  %
  % Inputs:
  %   P  #P by 2 list of positions (bbd should be ~40)
  %   E  #E by 3 list of edge indices into P
  % Outputs:
  %   V  #V by 3 list of positions of cutter part
  %   F  #F by 3 list of triangle indices on cutter part
  %
  % Example:
  % 
  %     clf;
  %     axis([-40 40 -40 40]);axis equal;
  %     P = get_pencil_curve();
  %     E = [1:size(P,1);2:size(P,1) 1]';
  %     [V,F] = cookie_cutter(P,E);
  %     tsurf(F,V);axis equal;
  %
  function [V,O] = robust_overlap(S,E)
    sqr_len = mean(sum((S(E(:,1),:)-S(E(:,2),:)).^2,2));
    max_area = sqr_len/9;
    [V,F] = cdt(S,E,'TriangleFlags',sprintf('-qYY -a%g',max_area),'Quiet',true);
    w = winding_number(S,E,barycenter(V,F));
    F = F(w>0.5,:);
    [V,I] = remove_unreferenced(V,F);
    F = I(F);
    O = outline(F);
    [V,I] = remove_unreferenced(V,O);
    O = I(O);
    [V,~,I] = remove_duplicate_vertices(V,1e-3);
    O = I(O);
    % remove degnerate edges
    O = O(O(:,1) ~= O(:,2),:);
    % remove duplice edges
    [~,I] = unique(sort(O,2),'rows');
    O = O(I,:);
    assert(size(unique(sort(O,2),'rows'),1) == size(O,1));
    [V,F] = triangle(V,O,[],'Flags','-YY','Quiet');
    w = abs(winding_number(V,O,barycenter(V,F)));
    F = F(w>0.5,:);
    O = outline(F);
    [V,I] = remove_unreferenced(V,O);
    O = I(O);
    assert(size(unique(sort(O,2),'rows'),1) == size(O,1));

  end

  % thin wall
  offset = 0.5;
  [SI,SO,EI,EO] = offset_curve(P,offset);
  S = [SI;SO];
  E = [EI;size(SI,1)+EO];
  [S,E] = robust_overlap(S,E);

  % thick wall
  k_offset = 2;
  [KI,KO,KEI,KEO] = offset_curve(P,k_offset);
  K = [KI;KO];
  KE = [KEI;size(KI,1)+KEO];
  [K,KE] = robust_overlap(K,KE);

  % height of cutting wall
  h = 15;
  % height of lip
  th = 3;
  % height of stamper
  sth = 5;
  [CV,CF] = wedding_cake({K,S},{KE,E},[th,h]);

  tol = 0.8;
  [sP,~,sE,~] = offset_curve(P,offset+tol);
  sE = fliplr(sE);
  [sP,sE] = robust_overlap(sP,sE);
  % center
  c = mean(sP);
  % center should be inside
  assert(abs(winding_number(sP,sE,c))>0.5,'Uhoh, center for post fell outside curve. Too non-convex');
  % radius determined by closest point on curve
  r = min(pdist2(K,c));
  r = 0.9*r;
  r = min(r,4);
  theta = linspace(0,2*pi,100)';
  theta = theta(1:end-1);
  hP = bsxfun(@plus,c,r*[cos(theta) sin(theta)]);
  hE = [1:size(hP,1);2:size(hP,1) 1]';
  [SV,SF] = wedding_cake({sP,hP},{sE,hE},[sth,2*h]);

  % concatenate
  F = [CF;size(CV,1)+SF];
  %V = [CV;bsxfun(@plus,[max(CV(:,1))*2 0 0],SV)];
  % Cutter needs to be flipped!
  V = [bsxfun(@times,[-1 1 1],CV);bsxfun(@plus,[abs(max(CV(:,1))-min(CV(:,1)))*1.1 0 0],SV)];
end
