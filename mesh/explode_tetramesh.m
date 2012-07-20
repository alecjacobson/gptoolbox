function h = explode_tetramesh(T,V,steps)
  % EXPLODE_TETRAMESH animated, interactive tetramesh plot of tets of a tet
  % mesh being pulled away from each 
  %
  % h = explode_tetramesh(T,V)
  %
  % Inputs:
  %   T  #T by 4 list of tetrahedra indices
  %   V  #V by dim by #{frames | 1} list of vertex positions
  %   steps  # of explosion animation steps {100}
  % Outputs:
  %   h  handle to plot
  %
  % See also: tsurf, trisurf, tetramesh, animated_trisurf, animated_tetramesh
  %

  % Remap vertices and tets so that all corners are separated
  V = [ ...
    V(T(:,1),:); ...
    V(T(:,2),:); ...
    V(T(:,3),:); ...
    V(T(:,4),:); ...
  ];
  T = [ ...
    0*size(T,1) + (1:size(T,1)); ...
    1*size(T,1) + (1:size(T,1)); ...
    2*size(T,1) + (1:size(T,1)); ...
    3*size(T,1) + (1:size(T,1)); ...
    ]';
  factor = 4;
  if ~exist('steps','var')
    steps = 100;
  end
  s = linspace(1,factor,steps);
  % centroid
  c = mean(V);
  VV = zeros(size(V,1),size(V,2),steps);
  for ii = 1:steps
    sii = s(ii);
    % dialate globally
    Vii = bsxfun(@plus,sii*V,-sii*c+c);
    % contract locally
    % centroids of each tet
    Cii = ( ...
      Vii(T(:,1),:) + ...
      Vii(T(:,2),:) + ...
      Vii(T(:,3),:) + ...
      Vii(T(:,4),:))/4;
    Vii = Vii/sii+repmat(Cii,4,1)*(1-1/sii);
    VV(:,:,ii) = Vii;
  end
  h = animated_tetramesh(T,VV);
  axis equal;
  A = reshape(axis,2,[]);
  % grow by factor about center
  dA = [A(1,:)-A(2,:);A(2,:)-A(1,:)];
  axis(reshape(A+dA,1,prod(size(A))));
  % don't resize mesh
  axis manual;
  
end
