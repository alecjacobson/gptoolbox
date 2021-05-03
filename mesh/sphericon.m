function [V,FV,U,FU,N,FN] = sphericon(r,n)
  % SPHERICON  Construct a mesh of a sphericon
  %
  % [V,F] = sphericon(r)
  %
  % Input:
  %   r  radius
  %   n  number of vertices along side
  % Outputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices
  %


  warning('alternative');

  U = [];
  y = 0;
  % minimum
  k = min(3,n);
  for j = 1:n
    y = (j-1)/(n-1);
    if j<k
      x = linspace(0,1,n)';
    else
      x = linspace(0,1,n-j+k)';
    end
    x(:,2) = y;
    U = [U;x];
  end
  sct(U,'.');
  F = delaunay(U.*[1 2]);
  F = F(sum(reshape(U(F,2),size(F))<1,2)>=2,:);
  A = U(:,1)*pi;
  R = 1-U(:,2);
  V = [ R.*[cos(A) sin(A)] U(:,2)];
  %[U,F] = cylinder_mesh(1,2*n,'Stacks',n);
  %F = F(floor( ((1:size(F,1))'-1)/(2*n)) < 2*n-1,:);
  %% get angle _before_ degnerat-izing
  %A = atan2(U(:,2),U(:,1));
  %V = [U(:,1:2).*(1-U(:,3)) U(:,3)];
  %F = F(barycenter(V(:,2),F)>0,:);
  %[V,~,J,F] = remove_unreferenced(V,F);
  %A = A(J);
  %U = U(J,:);
  N = [ sqrt(2)/2.*cos(A) sqrt(2)/2.*sin(A) repmat(sqrt(2)/2,size(A,1),1)];
  FV = F;
  FN = F;
  FU = F;
  FV = [FV;size(V,1)+FV];
  FN = [FN;size(N,1)+FN];
  FU = [FU;size(U,1)+FU];
  N = [N;-N(:,1) N(:,2) -N(:,3)];
  V = [V;-V(:,1) V(:,2) -V(:,3)];
  U = [U;U];
  FV = [FV;size(V,1)+FV];
  FN = [FN;size(N,1)+FN];
  FU = [FU;size(U,1)+FU];
  R = axisangle2matrix([0 0 1],pi)*axisangle2matrix([0 1 0],pi/2);
  N = [N;N*R];
  V = [V;V*R]*r;
  U = [U;U];
  % min(min(edge_lengths(V,FV))) == (pi*n^-2)
  [V,~,J] = remove_duplicate_vertices(V,(pi*n^-2)/10);
  FV = J(FV);
  return;

  % Half a bicone
  [CV,CF] = create_regular_grid(n,n,0,0);
  m = size(CF,1);
  % bottom left 
  BL = (2*(n-1)-mod((1:m)'-1, 2*(n-1)))<=(1+2*(floor(((1:m)'-1)/(2*(n-1)))));
  CF = flip_ears(CV,CF);
  % Rotate 45°
  R = [cos(pi/4) sin(pi/4);-sin(pi/4) cos(pi/4)];
  CV = bsxfun(@minus,CV*R*2/sqrt(2),[0,1]);

  % Not uniform
  %BV = [CV(:,1) real(sqrt((1-abs(CV(:,2))).^2-CV(:,1).^2)) CV(:,2)];
  % More uniform
  R = (1-abs(CV(:,2)));
  T = ((CV(:,1)+(1-abs(CV(:,2))))./((1-abs(CV(:,2)))*2)*2-1)*pi/2;
  T(isnan(T)) = 0;
  BV = [R.*sin(T) R.*cos(T) CV(:,2)];

  % Rotate and connect
  V = [BV;BV*axisangle2matrix([0 1 0],-pi/2)*axisangle2matrix([1 0 0],-pi)];
  % Rotate everything to lie flat and scale
  V = r*V*axisangle2matrix([0,1,0],pi/4);
  F = fliplr([CF;size(BV,1)+CF]);
  I = [BL+1;2+(BL+1)];
  [V,~,J] = remove_duplicate_vertices(V,r*0.1/n);
  F = J(F);

  %% Glue together
  %[V,F] = clean(V,F,'MinDist',1e-8,'MinArea',0,'MinAngle',0, ...
  %    'SelfIntersections','ignore','SmallTriangles','remove');
  %% double check that it's solid
  %statistics(V,F,'Fast',true)

end
