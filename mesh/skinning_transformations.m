function [varargout] = skinning_transformations(C,P,BE,new_C,RP)
  % SKINNING_TRANSFORMATIONS compute transformations for each handle of a
  % skinning deformation setup.
  %
  % [TR] = skinning_transformations(C,P,BE,new_C,RP)  
  % [T,R,S] = skinning_transformations(C,P,BE,new_C,RP) Same as above except
  %   rotations, scales and translations are separated, so that TR(:,1:dim,i) =
  %   R(:,:,i)*S(:,:,i) and TR(:,dim+1,:) = TR(:,3,:)= permute(T,[2,3,1])
  % [T,AX,AN,S,O] = skinning_transformations(C,P,BE,new_C,RP) Same as above
  %   except instead of representing rotations by a dim by dim matrix
  %   representing them as axes, AX, and angles, AN. For 2D control points C,
  %   AX(i,:) = [0 0 1]
  %
  % Inputs:
  %  C  #C by dim list of control vertex positions
  %  P  #P list of indices into C for point controls, { 1:size(C,1) }
  %  BE  #BE by 2 list of bones, pairs of indices into C, connecting control 
  %    vertices in C
  %  new_C  #C by dim list of new control point positions
  %  RP #dim by #dim by #P list of rotations applied to each control point in P
  %     or when dim == 2, can be #P list of angles
  % Outputs
  %  TR #dim by #dim+1 by #P+#BE list of transformations, where TR(1:2,3,i) is 
  %    the translation for control handle i and TR(1:2,1:2,i) is the 
  %    rotation/scale for control handle i
  %  or
  %  T #P+#BE by #dim list of translation vectors, where T(i,:) is the
  %    translation vector for control handle i
  %  R #dim by #dim by #P+#BE list of rotations matrices, where R(1:2,1:2,i) is
  %    the rotation matrix for control handle i
  %  S #dim by #dim by #P+#BE list of scale matrices, where S(1:2,1:2,i) is the
  %    scale matrix for control handle i
  %  or
  %  T #P+#BE by #dim list of translation vectors, where T(i,:) is the
  %    translation vector for control handle i
  %  AX #P+#BE by #dim list of rotation axis vectors, where AX(i,:) is the
  %    rotation axis vector for control handle i
  %  AN #P+#BE list of rotation angles, where AN(i) is the rotation angle for
  %    control handle i
  %  S #dim by #dim by #P+#BE list of scale matrices, where S(1:2,1:2,i) is the
  %    scale matrix for control handle i
  %  O #P+#BE + #dim list of origins
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: lbs, dualquatlbs

  % number o dimensions
  dim = size(C,2);

  if(exist('P','var'))
    if(size(P,1) == 1)
      P = P';
    end
  else
    P = (1:size(C,1))';
  end

  % point handles should be empty or a column vector
  assert(isempty(P) || (size(P,2) == 1));
  % BE should be empty or be #bones by 2 list of indices
  assert( isempty(BE) || (size(BE,2) == 2));
  % indices should be valid
  assert(isempty(P) || (max(P) <= size(C,1)));
  assert(isempty(BE) || max(BE(:)) <= size(C,1));

  % number of point controls
  np = numel(P);
  % number of bone controls
  nb = size(BE,1);


  assert(all(size(C) == size(new_C)));
  % if rotations are not given or given empty matrix just use identity rotation
  if(~exist('RP','var') || isempty(RP) || all(RP(:)==0))
    RP = repmat(eye(dim,dim),[1,1,np]);
    if(nargout == 5)
      ANp = zeros([np,1]);
      AXp = repmat([0 0 1],[np 1]);
    end
  elseif(numel(RP) == np)
    % user provided rotations as angles (for 2D only)
    assert(dim == 2);
    % copy angles to new name
    ANp = reshape(RP,np,1);
    AXp = repmat([0 0 1],[np 1]);
    % convert to matrices
    RP = zeros([dim dim np]);
    RP(1,1,:) =  cos(ANp);
    RP(1,2,:) = -sin(ANp);
    RP(2,1,:) =  sin(ANp);
    RP(2,2,:) =  cos(ANp);
  else
    if(nargout == 5)
      assert(false);
    end
  end
  % should have rotation for each control point
  assert(size(RP,1) == dim);
  assert(size(RP,2) == dim);
  assert(size(RP,3) == np);

  % allocate space for transformations
  R = repmat(eye(dim,dim),[1,1,np+nb]);
  S = repmat(eye(dim,dim),[1,1,np+nb]);

  % insert point rotations (and scales)
  R(1:dim,1:dim,1:np) = RP;
  % for now, don't allow scaling at points
  %S(1:dim,1:dim,1:np) = RS;

  % insert bone rotations (and scales)
  if(nb > 0)
    % edge vector at rest state
    rest = C(BE(:,2),:) - C(BE(:,1),:);
    rest_norms = normrow(rest);
    % edge vector at pose
    pose = new_C(BE(:,2),:) - new_C(BE(:,1),:);
    pose_norms = normrow(pose);

    % first take care of scaling anisotropically along bone
    % rotate so that rest vector is x-axis
    XR = zeros([dim dim nb]);
    % signed angle that takes rest to x-axis
    phi = -atan2(rest(:,2),rest(:,1));
    % rotation is this angle around axis of cross between these two vectors
    if(dim == 2)
      % in 2d rotations are always around "z-axis" and rotation matrix is
      % simple
      XR(1,1,:) =  cos(phi);
      XR(1,2,:) = -sin(phi);
      XR(2,1,:) =  sin(phi);
      XR(2,2,:) =  cos(phi);
    else
      %axis = cross(repmat([1 0 0],2),rest);
      phi_axis = [0*rest(:,1) -rest(:,3) rest(:,2)];
      % convert axis,angle to rotation matrix
      % 3D not supported yet
      assert(false);
    end

    % store rotation in scale spot
    S(:,:,(np+1):(np+nb)) = XR;
    % scale along x-axis
    S(1,:,(np+1):(np+nb)) = ...
      XR(1,:,:) .* ...
      permute(repmat((pose_norms./rest_norms),[1 dim 1]),[3 2 1]);
    % unrotate each scale
    for( ii = 1:nb)
      S(:,:,np+ii) = XR(:,:,ii) \ S(:,:,np+ii);
    end


    % rotation is this angle around axis of cross between these two vectors
    if(dim == 2)
      % in 2d rotations are always around "z-axis" and rotation matrix is
      % simple

      % use 2D signed angle
      ANb = atan2(pose(:,2),pose(:,1)) - atan2(rest(:,2),rest(:,1));
      R(1,1,(np+1):(np+nb)) =  cos(ANb);
      R(1,2,(np+1):(np+nb)) = -sin(ANb);
      R(2,1,(np+1):(np+nb)) =  sin(ANb);
      R(2,2,(np+1):(np+nb)) =  cos(ANb);
      AXb = repmat([0 0 1],[nb 1]);
    else
      % angle between rest and pose
      ANb = acos(dot(rest,pose,2) ./ (rest_norms .* pose_norms));
      % treat degenerate angles as 0
      ANb(isnan(ANb)) = 0;
      AXb = cross(rest,pose,2);
      % convert axis,angle to rotation matrix
      assert(false);
    end

  else
    ANb = [];
    AXb = [];
  end

  % add to translation, differnce in origin and rotation/scale applied to 
  % apply scale before rotation
  RS = zeros(dim*(np+nb),dim);
  for( ii = 1:(np+nb))
    RS((1+dim*(ii-1)):(dim*ii),1:dim) = R(1:dim,1:dim,ii) * S(:,:,ii);
  end
  % origin
  T = zeros(np+nb,dim);
  % origins
  O = zeros(np+nb,dim);
  % origins for point is simply the point at rest state
  O(1:np,:) = C(P,:);
  % insert point translations
  T(1:np,:) = new_C(P,:) - C(P,:);
  if(nb > 0)
    % origin for bones, just use start point
    O((np+1):(np+nb),:) = C(BE(:,1),:);
    % insert bone translations, just use start point
    T((np+1):(np+nb),:) = new_C(BE(:,1),:) - C(BE(:,1),:);
  end

  % assemble output based on number of output parameters
  if(nargout <= 1)
    % linear transformations
    T = (T+O-stacktimes(RS,O));
    TR = zeros([dim dim+1 (np+nb)]);
    TR(:,1:dim,:) = permute(reshape(RS',[dim dim (np+nb)]),[2 1 3]);
    TR(:,dim+1,:)= permute(T,[2,3,1]);
    varargout{1} = TR;
  elseif(nargout == 3)
    % translations and rotations and scales
    T = (T+O-stacktimes(RS,O));
    varargout{1} = T;
    varargout{2} = R;
    varargout{3} = S;
  elseif(nargout == 5)
    %T = (T+O-stacktimes(R,O));
    % translations axis angle scale origin
    R = reshape(permute(R,[2 1 3]),[dim dim*(np+nb)])';
    T = (T+O-stacktimes(R,O));
    varargout{1} = T;
    varargout{2} = [AXp;AXb];
    varargout{3} = [ANp;ANb];
    varargout{4} = S;
    varargout{5} = O;
  end

end
