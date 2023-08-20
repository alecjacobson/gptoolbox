function [phi,F,I,C,fold_fun] = orbifold(Ubar,Fbar,Cbar,orbi_type)
  % ORBIFOLD Compute a conformal, globally injective parameterization of a
  % sphere-topology surface, given selected cone singularities and type.
  % "Orbifold Tutte Embedding" Aigerman and Lipman.
  %
  % [phi,F,I,C] = orbifold(Ubar,Fbar,Cbar,orbi_type)
  %
  % Inputs:
  %   Ubar  #Ubar by 3 list of input sphere-topology mesh vertex positions
  %   Fbar  #Fbar by 3 list of triangle indices into Ubar
  %   Cbar  3-long list of cone singularity indices into rows of Ubar {found
  %     using farthest_points}
  %   orbi_type  followed by orbifold type {'i'}
  % Outputs:
  %   phi  #phi by 2 list of 2D UV positions
  %   F  #F by 3 list of triangle indices into rows of phi
  %   I  #phi list of indices so that U = Ubar(I,:)
  %   C  4-long list of cone singularity indices into Cbar (with repeated
  %     indices)
  %   fold_fun  function so that fold_fun(X) maps points outside [0,1]² to
  %     appropriate tiling-modulo position in [0,1]²
  %

  if ~exist('Cbar','var') || isempty(Cbar)
    [~,Cbar] = farthest_points(Ubar,3);
    Cbar = Cbar([1 3 2]);
  end
  if ~exist('orbi_type','var') || isempty(orbi_type)
    orbi_type = 'i';
  end
  assert(strcmp(orbi_type,'i'));

  % Determine cut edges. This don't matter but don't let them cross or get too
  % close.
  Abar = adjacency_matrix(Fbar);
  Gbar = graph(Abar);
  P12 = shortestpath(Gbar,Cbar(1),Cbar(2));
  Ebar = [P12(1:end-1);P12(2:end)]';
  block = (Abar * accumarray(P12(2:end-1)',1,[size(Ubar,1) 1]))>0;
  block(Cbar) = 0;
  Abar(block,:) = 0;
  Abar(:,block) = 0;
  Gbar = graph(Abar);
  P23 = shortestpath(Gbar,Cbar(2),Cbar(3));
  Ebar = [Ebar;[P23(1:end-1);P23(2:end)]'];
  % Which edge (if any) is each original vertex on?
  J = zeros(size(Ubar,1),1);
  J(P12(2:end-1)) = 1;
  J(P23(2:end-1)) = 2;

  % Cut the mesh along these edges
  [F,I] = cut_edges(Fbar,Ebar);
  U = Ubar(I,:);
  C = [find(I==Cbar(1));find(I==Cbar(2));find(I==Cbar(3))];


  % Recover boundary paths. Rather than try to sort things out with I just
  % retrace shortest paths along the (only new) boundary, this ensures we choose
  % the same side through singularities
  dVmC = find(J(I) > 0);
  dV = [dVmC;C];
  VmdV = setdiff((1:size(U,1))',dV);
  A = adjacency_matrix(F);
  A(VmdV,:) = 0;
  A(:,VmdV) = 0;
  G = graph(A);
  % This is specific to types i,ii,iii
  P12_left = shortestpath(G,C(1),C(2));
  P23_left = shortestpath(G,C(2),C(4));
  P12_right = shortestpath(G,C(1),C(3));
  P23_right = shortestpath(G,C(3),C(4));
  IJ = { 
    [P12_left(2:end-1); P12_right(2:end-1)]', 
    [P23_left(2:end-1); P23_right(2:end-1)]'};

  L = cotmatrix(U,F);
  % Could also use intrinsic. This might matter when strict safety is needed.
  % This is expensive so it should be an option.
  %L = intrinsic_delaunay_cotmatrix(U,F);

  % Build linear system

  % First, for all interior vertices vi ∈ V \ ∂V we set the standard discrete harmonic equation as in (2),
  AVmdV = repdiag(L(VmdV,:),2);
  AC = repdiag(sparse(1:numel(C),C,1,numel(C),size(U,1)),2);

  orbi_type = 'i';
  switch orbi_type
  case 'i'
    % Cbar 1→2→3
    % C 1→2→4 and 1→3→4
    % consider each boundary segment pair
    R = cat(3,[0 -1;1 0],[0 1;-1 0]);
    r = [pi/2 pi pi pi/2];
    XC = [0 0;0 1;1 0;1 1];
    % Which vertices are we rotating about?
    CI = [C(1) C(4)];
  end


  AdVmC = [];
  for k = 1:size(R,3)
    Rk = R(:,:,k);
    Ik = IJ{k}(:,1);
    Jk = IJ{k}(:,2);
    AIk = repdiag(L(Ik,:),2) + kron(Rk,L(Jk,:));
    nk = numel(Ik);
    DC = @(Ik) sparse( ...
      [1:nk;1:nk]', ...
      [Ik repmat(CI(k),nk,1)], ...
      repmat([1 -1],nk,1), ...
      nk,size(U,1));
    AJk = kron(Rk,DC(Jk)) - repdiag(DC(Ik),2);
    AdVmC = [AdVmC;AIk;AJk];
  end

  A = [AVmdV;AdVmC;AC];
  b = [zeros(2*(numel(VmdV)+numel(dVmC)),1);XC(:)];
  % is A(VmC,VmC) secretly symmetric?
  phi = reshape(A\b,size(U,1),2);

  fold_fun = @(phi0) fold(phi0,r,phi(C,:));

  function phi = fold(phi0,r,pC)
    phi = phi0;
  
    Rot = @(r) [cos(r) -sin(r);sin(r) cos(r)];
    left = phi(:,1)<0 & (phi(:,2)>=0 & phi(:,2)<=1);
    right = phi(:,1)>1 & (phi(:,2)>=0 & phi(:,2)<=1);
    top = phi(:,2)>1 & (phi(:,1)>=0 & phi(:,1)<=1);
    bottom = phi(:,2)<0 & (phi(:,1)>=0 & phi(:,1)<=1);
    top_left = phi(:,1)<0 & phi(:,2)>1;
    bottom_right = phi(:,1)>1 & phi(:,2)<0;
  
    phi(left,:)   = phi(left,:)*Rot(r(1));
    phi(bottom,:) = phi(bottom,:)*Rot(r(1))';
    phi(right,:)  = (phi(right,:)-pC(4,:))*Rot(r(4))+pC(4,:);
    phi(top,:)    = (phi(top,:)-pC(4,:))*Rot(r(4))' +pC(4,:);
    phi(top_left,:)   = (phi(top_left,:)-pC(2,:))*Rot(r(2))+pC(2,:);
    phi(bottom_right,:)   = (phi(bottom_right,:)-pC(3,:))*Rot(r(3))+pC(3,:);
  
    %clf;
    %hold on;
    %tsurf([1 2;2 3;3 4;4 1],[0 0;1 0;1 1;0 1]);
    %sct(phi);
    %hold off;
    %error
  end
end
