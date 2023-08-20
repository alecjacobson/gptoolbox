function tsh = plot_taper(V,E,R,varargin)
 % tsh = plot_taper(V,E,R,varargin)
 % 
 % Example:
 %
 % % (P,C) cubic spline with radii R at control points
 % [tV,tE] = spline_to_poly([P R],C,0.1);
 % plot_taper(tV(:,1:2),tE,tV(:,3),'FaceColor','r','EdgeColor','none');



  if numel(R) == 1
    R = repmat(R,size(V,1),1);
  end
  % Circles at each vertex
  th = linspace(0,2*pi,50);
  th = th(1:end-1)';
  CV = [cos(th) sin(th)];
  CE = [1:size(CV,1); 2:size(CV,1) 1]';
  [CV,CF] = triangulate(CV,CE);
  assert(size(V,1) == size(R,1) || numel(R) == 1);
  J = unique(E(:));
  [CV,CF,~,I] = repmesh(CV,CF,V(J,:),R(J));

  % parallelograms at each edge
  EV = V(E(:,2),:)-V(E(:,1),:);
  N = normalizerow(EV*[0 -1;1 0]);
  PV = [ ...
    V(E(:,1),:) + R(E(:,1)).*N ; ...
    V(E(:,1),:) + -R(E(:,1)).*N ; ...
    V(E(:,2),:) + -R(E(:,2)).*N ; ...
    V(E(:,2),:) + R(E(:,2)).*N];
  PQ = (1:size(E,1))' + (0:3)*size(E,1);
  PF = [PQ(:,1:3);PQ(:,[1 3 4])];

  VV = [CV;PV];
  FF = [CF;size(CV,1)+PF];
  tsh = tsurf(FF,VV,varargin{:});


end
