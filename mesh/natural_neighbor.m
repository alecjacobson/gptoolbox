function [W] = natural_neighbor(V,C)
  % NATURAL_NEIGHBOR  computes Sibson's AKA natural neighbor coordinates at
  % points in V for controls (samples) C
  %
  % [W] = natural_neighbor(V,C)
  %
  % Inputs:
  %  V  #V by 2 list of domain positions
  %  C  #C by 2 list of control positions
  %
  % Output:
  %  W  #V by #C list of natural neighbor coordinates
  %

  % number of domain points
  n = size(V,1);
  % number of original samples
  m = size(C,1);
  assert(size(V,2) == 2);
  assert(size(V,2) == size(C,2));
  
  % large number
  E = 1e10*max(abs(V(:)));
  % enclose samples with ghost points
  %G = [-E -E; -E E; E -E; E E; E 0;0 E;-E 0; 0 -E];
  ng = 10*m;
  theta = linspace(0,2*pi,ng+1)';
  theta = theta(1:(end-1));
  G = E.*[cos(theta) sin(theta)];

  % append ghost points to samples
  C = [C;G];

  % create delaunay triangulation of samples and ghost points
  DT = DelaunayTri(C);
  triplot(DT)

  W = zeros(n,size(C,1));
  % compute natural neighbor interpolation for each control point
  for ii = 1:size(C,1)
    progressbar(ii,size(C,1));
    % function value
    b = zeros(size(C,1),1);
    b(ii) = 1;
    % create interpolation table
    Fnat = TriScatteredInterp(DT,b,'natural');
    W(:,ii) = Fnat(V(:,1),V(:,2));
  end
  % only keep coordinate for true samples
  W = W(:,1:m);
  % renormalize
  W = W./repmat(sum(W,2),1,m);
end
