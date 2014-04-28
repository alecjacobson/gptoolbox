function [U] = stiff_points_lbs(V,C,P,SE,new_C,WP,WS)
  % STIFF_POINTS_LBS Linear blend skinning with stiff points
  %
  % [U] = stiff_points_lbs(V,C,P,BE,new_C,WP,WS)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   C  #C by dim list of control vertex rest positions
  %   P  #P list of indice into C of point controls
  %   SE  #B by 2 list of stiffening bone indices into *P*
  %   new_C #C by dim list of update postions for control vertices
  %   W  #V by #P list of point weights
  %   WS  #V by #S list of bone endpoint weights cooresponding to BE(:,2)
  % Output:
  %   U  #V by dim list of new vertex positions
  % 


  assert(size(V,2) == size(C,2));
  % number of dimensions
  dim = size(V,2);
  % number of domain vertices
  n = size(V,1);
  % number of points
  np = numel(P);
  % number of stiffening edges
  ns = size(SE,1);
  assert(max(P) <= size(C,1));
  assert(min(P) >= 1);
  assert(size(SE,2) == 2);
  assert(max(SE(:)) <= np);
  assert(min(SE(:)) >= 1);
  assert(size(WP,1) == n);
  assert(size(WP,2) == np);
  assert(size(WS,1) == n);
  assert(size(WS,2) == ns);
  assert(all(size(C) == size(new_C)));

  % 
  % u = ∑_i wpi (v) [ pi' + Ri(v) * ( v - pi ) ]
  % where
  % Ri(v) = ∑_j∈N(i) wsj(v) * Rij
  % and Rij is the bone rotation of points i and j

  % Mask weights with bone incidence
  % M(i,k) is 1 iff handle i is incident on bone k
  M = sparse([SE(:,1)' SE(:,2)'],[1:ns 1:ns],1);
  % M(j,i,k) is 1 only iff handle i is incident on bone k
  M = permute(repmat(full(M),[1 1 n]),[3 1 2]);
  WSP = permute(repmat(WS,[1 1 np]),[1 3 2]) .* M;
  % normalize per vertex per handle
  WSP = WSP ./ repmat(sum(WSP,3),[1 1 ns]);
  % Why do I have to do this?
  WSP(isnan(WSP)) = 0;

  rest = C(P(SE(:,2)),:) - C(P(SE(:,1)),:);
  pose = new_C(P(SE(:,2)),:) - new_C(P(SE(:,1)),:);

  % compute bone rotations
  [w,a] = axisanglebetween(pose,rest,[ 0 0 1]);

  if(dim == 2)
    % convert to signed angle
    a = a.*w(:,3);
    % A(i,k) rotation of bone k seen by handle i
    A = sparse([SE(:,1)' SE(:,2)'],[1:ns 1:ns],[a a]);
    % A(j,i,k) rotation of bone k seen by handle i weighted according to vertex
    % j of domain
    A = permute(repmat(full(A),[1 1 n]),[3 1 2]);
    %WSP = permute(repmat(WS,[1 1 np]),[1 3 2]);

    % perform weighted average
    A = sum(WSP .* A,3);

    % convert to matrix
    R = zeros([dim dim n np]);
    R(1,1,:,:) = cos(A);
    R(1,2,:,:) = -sin(A);
    R(2,1,:,:) = sin(A);
    R(2,2,:,:) = cos(A);
  else
    % convert to quaternions
    Q = axisangle2quat(w,a);
    assert(false);
  end
  % should now have a rotation for each domain vertex for each handle

  % repeat weights for each dimension
  WPdim = permute(repmat(WP,[1 1 dim]),[1 3 2]);
  U = zeros([n dim np]);
  for ii = 1:np
    Rii = reshape(R(:,:,:,ii),[ dim dim*n])';
    U(:,:,ii) = WPdim(:,:,ii).*( ...
      repmat(new_C(P(ii),:),[n 1]) + ...
      stacktimes(Rii,V - repmat(C(P(ii),:),[n 1])));
  end
  U = sum(U,3);




  %% subtract handle rest position from each domain vertex, VmC(j,:,i) gives
  %% V(j,:) - C(i,:)
  %VmC = repmat(V,[1 1 np]) - permute(repmat(C(P,:),[1 1 n]),[3 2 1]);

  %R = reshape(permute(R,[3 1 2 4]),[dim*n,dim,np]);
  %R
  %RVmC = zeros([n dim np]);
  %% RVmC(j,:,i) = R(:,:,j,i) * VmC(j,:,i)
  %for ii = 1:np
  %  RVmC(:,:,ii) = stacktimes(R(:,:,ii),VmC(:,:,ii));
  %end

  %% repeat weights for each dimension
  %WPdim = permute(repmat(WP,[1 1 dim]),[1 3 2]);
  %% repeat displacement for each vertex
  %new_CV = permute(repmat(new_C,[1 1 n]),[3 1 2]);
  %% weighted average
  %U = sum(WPdim.*(new_CV + RVmC),3);
end
