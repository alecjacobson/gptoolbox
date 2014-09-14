function [W,PP] = bone_heat(V,F,C,P,BE,CE)
  % BONE_HEAT computes "bone heat" weights from "Automatic Rigging and
  % Animation of 3D Characters" by Baran and Popovic.
  %
  % W = bone_heat(V,F,C,P,BE,CE)
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices, for 3D F is #F by 4, for 2D F is #F by 3
  %  C  list of control vertex positions
  %  P  list of indices into C for point controls, { 1:size(C,1) }
  %  BE  list of bones, pairs of indices into C, connecting control vertices, 
  %    { [] }
  %  CE  list of "cage edges", pairs of indices into ***P***, connecting
  %    control ***points***. A "cage edge" just tells point boundary conditions 
  %    to vary linearly along straight lines between the end points, and to be
  %    zero for all other handles. { [] }
  %
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  %

  assert(isempty(P) || (prod(size(P)) == max(size(P))));
  % E should be empty or be #bones by 2 list of indices
  assert( isempty(BE) || (size(BE,2) == 2));
  assert( isempty(CE) || (size(CE,2) == 2));

  nce = size(CE,1);
  % number of mesh vertices
  n = size(V, 1);
  % number of point controls
  np = numel(P);
  % number of bone edges
  nb = size(BE,1);
  m = np+nb;

  % determine if we can use accelerated visibility test
  if 3==exist('bone_visible_embree','file')
    bone_visible_func = @bone_visible_embree;
  else
    bone_visible_func = @bone_visible;
    warning('Using non accelerated visibility test');
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % DETECTING VISIBILITY IS BY FAR THE BOTTLENECK: 99% of time is spent in
  % point_visible and bone_visible
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % build visibility mask, vis_mask, so that vis_mask(j,i) says whether shortest
  % line segment from handle i to mesh vertex j is inside mesh
  vis_mask = true(n,m);
  % build cage edge visibility mask, same as above but for the cage edges
  CE_vis_mask = true(n,nce);

  tic;
  % build distance vector, D, so that D(j,i) is the length of the shortest line
  % segment from handle i to mesh vertex j
  D = zeros(n,m);
  % collect point distances and visibility
  for ii = 1:np
    fprintf('Point %d out of %d\n',ii,np);
    D(:,ii) = sqrt(sum((V - repmat(C(P(ii),:),[n 1])).^2,2));
    vis_mask(:,ii) = point_visible(V,F,C(P(ii),:));
    toc;
  end

  TCE = zeros(n,nce);
  DCE = zeros(n,nce);
  % collect cage edge distances and visibility
  for ii = 1:nce
    fprintf('Cage edge %d out of %d\n',ii,nce);
    % source point
    ps = C(P(CE(ii,1)),:);
    % dest point
    pd = C(P(CE(ii,2)),:);
    % project to lines 
    [t,sqr_d] = project_to_lines(V,ps,pd);
    % snap to end points
    sqr_d(t<0) = 1e16;
      %sum((V(t<0,:) - repmat(C(P(CE(ii,1)),:),[size(V(t<0,:),1) 1])).^2,2);
    sqr_d(t>1) = 1e16;
      %sum((V(t>1,:) - repmat(C(P(CE(ii,1)),:),[size(V(t>1,:),1) 1])).^2,2);
    TCE(:,ii) = t;
    DCE(:,ii) = sqrt(sqr_d);
    CE_vis_mask(:,ii) = bone_visible_func(V,F,ps,pd);
    toc;
  end

  % collect bone distances and visibility
  for ii = 1:nb
    fprintf('Bone %d out of %d\n',ii,nb);
    % source point
    ps = C(BE(ii,1),:);
    % dest point
    pd = C(BE(ii,2),:);
    % project to lines 
    [t,sqr_d] = project_to_lines(V,ps,pd);
    % snap to end points
    sqr_d(t<0) = ...
      sum((V(t<0,:) - repmat(C(BE(ii,1),:),[size(V(t<0,:),1) 1])).^2,2);
    sqr_d(t>1) = ...
      sum((V(t>1,:) - repmat(C(BE(ii,2),:),[size(V(t>1,:),1) 1])).^2,2);
    D(:,ii+np) = sqrt(sqr_d);
    vis_mask(:,ii) = bone_visible_func(V,F,ps,pd);
    toc;
  end

  % build matrix PP, so that PP(j,i) is 1 if D(j,i) < D(j,k) any k ??? i
  % unless cage edges exist then if D(j,ce(i)) is smallest then PP(j,ce(i,1)
  % and PP(j,ce(i,2) get 1-t and t respectively, t being the projection onto
  % the cage edge vector
  [minD,closest] = min([D DCE],[],2);
  % if cage edge was closest then we'll deal with that later
  indices = 1:n;
  % separate indices into those whose closest "handle" is a point or bone and
  % those whose closest "handles" is a cage edge
  indices_PB = indices(closest<=m);
  indices_CE = indices(closest>m);

  % list closest points/bones corresponding to mesh vertices whose closest
  % "handle" is a point/bone
  closest_PB = closest(closest<=m);
  % list closest cage edges corresponding to mesh vertices whose closest
  % "handle" is a cage edge
  closest_CE = closest(closest>m) - m;

  % point/bone contribution to PP
  PP_PB = sparse(indices_PB,closest_PB,1,n,m);

  if(isempty(CE))
    % no cages then only contribution is from points/bones
    PP = PP_PB;
  else
    % For points whose closest "handle" is a cage edge find their projection
    % onto that cage edge
    t = TCE(sub2ind(size(TCE),indices_CE,closest_CE'));
    % Cage edge contribution to PP
    PP_CE = sparse( ...
      [indices_CE indices_CE], ...
      [CE(closest_CE,1)' CE(closest_CE,2)'], ...
      [1-t t], ...
      n, ...
      m);
    % combine with points/bone contribution
    PP = PP_PB + PP_CE;
  end

  % build diagonal H matrix so that H(j,j) is 1/min(D(j,:),2))^2 if vis_mask(j,i)
  % otherwise 0
  closest_visible = true(n,1);
  % points/bones visibility, 1 if closest "handle" (point/bone) is visible
  closest_visible(indices_PB) = ...
    vis_mask(sub2ind(size(vis_mask),(indices_PB)',closest_PB));
  if(~isempty(CE))
    % cage edge visibility, 1 if closest "handle" (cage edge) is visible
    closest_visible(indices_CE) = ...
      CE_vis_mask(sub2ind(size(CE_vis_mask),(indices_CE)',closest_CE));
  end

  %only_boundary = false;
  %if dim == 2 && only_boundary
  %  % only use visibility info along outline
  %end

  H = sparse( ...
    1:n, ...
    1:n, ...
    (minD.^-2) .* closest_visible);
  % max-out min distances to avoid infs
  H = diag(min(diag(H),1e10));

  % build laplacian matrix 
  L = cotmatrix(V,F);
  % build mass matrix
  M = massmatrix(V,F,'voronoi');
  % make sure all rows are non-zero
  assert(all(sum(abs(M),2)~=0))
  % discrete laplace beltrami operator
  K = 0.5*(L);

  %% prefactor (-K+H) using lu decomposition, should be able to use cholesky...
  %[pfL,pfU,pfP,pfQ,pfR] = lu(-K+H);

  %W = zeros(n,m);
  %% loop over handles
  %for ii = 1:m
  %  % solve (-L+H) wi = H * PP(:,i)
  %  rhs = H * PP(:,ii);
  %  %W(:,ii) = pfQ*(pfU\(pfL\(pfP*(pfR\rhs)))); 
  %  W(:,ii) = (-K+H)\rhs;
  %end
  W = (-K+M*H)\(M*H*full(PP));

end
