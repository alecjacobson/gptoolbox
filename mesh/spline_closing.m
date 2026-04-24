function [cP,cC,cF] = spline_closing(P,C,sigma,nu)

  function plot_colored_splines(P,C,I)
    CM = cbrewer('Set1',10);
    for i = 1:size(C,1)
      c = mod(I(i)-1,size(CM,1))+1;
      plot_cubic(P(C(i,:),:),[],[],{{'MarkerFaceColor',CM(c,:)},{'Color',CM(c,:)}});
    end
  end

  function [C,sigma,M,N] = circle_center(A,B,sigma)
  % Compute center of circle C with radius sigma that goes through A and B
  % Triangle ABC has clockwise orientation.
  %
  % Inputs:
  %   A      #N by 2 list of A points
  %   B      #N by 2 list of B points
  %   sigma  scalar or #N by 1 list of radii
  %
  % Outputs:
  %   C      #N by 2 list of circle centers
  %   sigma  adjusted radii (clamped to >= |AB|/2)

  AB = B - A;
  d = sqrt(sum(AB.^2,2));

  % Ensure radius large enough
  sigma = max(sigma, d/2);

  % Midpoint
  M = (A + B)/2;

  % Distance from midpoint to center
  h = sqrt(max(sigma.^2 - (d/2).^2,0));

  % Unit perpendicular to AB
  % For AB = (x,y), clockwise normal is (y,-x)
  N = [AB(:,2), -AB(:,1)];
  N = N ./ d;

  % Center
  C = M + h .* N;
end


  minP = min(P(C,:));
  maxP = max(P(C,:));
  bbd = normrow(maxP-minP);

  [U,I,t,E] = spline_uniformly_sample(P,C,nu);
  [~,area] = centroid(U,E);
  if area<0
    P = P.*[1 -1];
    [cP,cC,cF] = spline_closing(P,C,sigma,nu);
    cP = cP.*[1 -1];
    return;
  end
  [~,area] = centroid(U,E);
  assert(area>0);

  N = normalizerow(cell2mat(arrayfun(@(i) cubic_tangent(P(C(I(i),:),:),t(i)),(1:size(I,1))','UniformOutput',false))*[0 -1;1 0]);
  Us = U + sigma*N;
  [~,Is,ts,K] = point_spline_squared_distance(Us,P,C);
  closed = normrow(U-K)>bbd*1e-8;

  closed_dst = closed(E(:,1)) & ~closed(E(:,2));
  closed_src = ~closed(E(:,1)) & closed(E(:,2));
  assert(~any(closed_dst & closed_src));
  assert(sum(closed_dst) == sum(closed_src));

  % transport all closed dst to src
  Xsrc = U(E(closed_src,2),:);
  Xdst_src = K(E(closed_dst,1),:);

  Xdst = U(E(closed_dst,1),:);
  Xsrc_dst = K(E(closed_src,2),:);
  D = pdist2(Xsrc,Xdst_src) + pdist2(Xsrc_dst,Xdst);
  [M,unsrc,undst] = matchpairs(D,1e26);
  assert(isempty(undst) && isempty(unsrc));
  % For each of the edges in closed_dst and closed_src
  %E_src = E(find(closed_src),:);
  %E_dst = fliplr(E(find(closed_dst),:));
  %E_src_dst = [E_src;E_dst];
  %[I_src,t_src] = find_closing_transition(P,C,sigma,E_src_dst,I,t,U,Is,ts,K);


  I_src = I(E(closed_src,2));
  t_src = t(E(closed_src,2));
  I_dst = I(E(closed_dst,1));
  t_dst = t(E(closed_dst,1));
  [sP,sC,sI,J] = spline_split(P,C,[I_src;I_dst],[t_src;t_dst]);
  iA = sC(J(1:end/2),4);
  iB = sC(J(end/2+1:end),4);
  J(1:end/2) = J(1:end/2) + 1;

  A = Xsrc(M(:,1),:);
  B = Xdst(M(:,2),:);
  AB = B-A;
  Q = circle_center(A,B,sigma);
  cP = [];
  cC = [];
  for i = 1:size(A,1)
    if any(isnan(Q(i,:)),'all')
      [cPi,cCi] = poly_to_spline([A(i,:);B(i,:)],[1 2]);
    else
      [cPi,cCi] = arc_to_cubics(A(i,:),B(i,:),Q(i,:));
    end
    assert(~any(isnan(cPi),'all'));
    cCi = cCi + size(sP,1) + size(cP,2);
    cCi(1,1) = iA(M(i,1));
    cCi(end,end) = iB(M(i,2));
    cC = [cC cCi'];
    cP = [cP cPi'];
  end
  cP = [sP;cP'];
  cC = cC';

  th = linspace(0,2*pi,100)';
  cV = [cos(th),sin(th)];
  cE = [1:numel(th); 2:numel(th),1]';
  [cV,cE] = repmesh(cV,cE,Q,sigma);

  sC = sC(~ismember(1:size(sC,1),J),:);
  cF = [true(size(sC,1),1); false(size(cC,1),1)];
  cC = [sC;cC];
  cc = connected_components(cC(:,[1 4]));
  val = accumarray(reshape(cC(:,[1 4]),[],1),1,[size(cP,1) 1]);
  bad = unique(cc(val==1));
  keep = ~any(ismember(cc(cC(:,[1 4])),bad),2);
  cC = cC(keep,:);
  cF = cF(keep);
  [cP,~,~,cC] = remove_unreferenced(cP,cC);

  %%plot_spline(P,C);
  %clf;
  %hold on;
  %%plot_spline(sP,sC);
  %sct(cP(val==1,:),'filled','b');
  %%plot_colored_splines(P,C,(1:size(C,1))');
  %%off = normrow(max(P)-min(P))*[1 0];
  %%plot_colored_splines(sP+off,sC,sI);
  %%sct(U(E(closed_src,2),:)+off,'g','filled');
  %%sct(Xdst,'r','filled','SizeData',100);
  %%sct(sP(iA,:),'filled','b');
  %%sct(sP(iB,:),'filled','c');

  %%sct(Xdst,'g','filled');
  %%tsurf(E_src_dst,U,'LineWidth',4,'EdgeColor','b');

  %%qvr(Xdst,Xdst_src-Xdst,0,'g');
  %%qvr(Xsrc,Xsrc_dst-Xsrc,0,'r');
  %%qvr(Xsrc(M(:,1),:),Xdst(M(:,2),:)-Xsrc(M(:,1),:),0,'b');
  %%sct(Q,'filled','b');
  %%tsurf(cE,cV);
  %plot_spline(cP,cC);

  %%txt(Xdst_src,[],'BackgroundColor',[0.7 0.99 0.7]);
  %%txt(Xsrc,[],'BackgroundColor',[0.99 0.7 0.7]);
  %%txt(Xdst,[],'BackgroundColor',[0.7 0.99 0.7]);
  %%txt(Xsrc_dst,[],'BackgroundColor',[0.99 0.7 0.7]);
  %hold off;
  %axis equal;

end
