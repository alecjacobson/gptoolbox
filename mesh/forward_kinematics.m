function [Q,T] = forward_kinematics(C,BE,P,dQ)
  % Compute absolute transformations from FK hierarchy and relative bone
  % rotations.
  %
  % Inputs:
  %   C  #C by dim list of joint positions
  %   BE  #BE by 2 list of bone edge indices
  %   P  #BE list of parent indices into BE
  %   dQ  #BE list of relative rotations
  % Outputs:
  %   Q  #BE list of absolute rotations
  %   T  #BE list of absolute translations

  function fk_helper(b)
    %function w = quatrotate(q,v)
    %  vv = [v 0];
    %  w = quatmultiply(quatmultiply(q,v),quatconj(q));
    %  w = w(1:3);
    %end
    if ~computed(b)
      p = P(b);
      if p < 1
        % base case for roots
        Q(b,:) = dQ(b,:);
        r = C(BE(b,1),:);
        T(b,:) = r - quatrotate(dQ(b,:),r);
      else
        % otherwise first compute parent's
        fk_helper(p);
        Q(b,:) = quatmultiply(Q(p,:),dQ(b,:));
        r = C(BE(b,1),:);
        T(b,:) = T(p,:) - quatrotate(Q(b,:),r) + ...
          quatrotate(Q(p,:),r);
      end
      computed(b) = true;
    end
  end


  % number of bones
  m = size(BE,1);
  assert(m == size(dQ,1));
  assert(m == size(P,1));

  % Whether we've already computed FK for this bone: dynamic programming
  computed = false(m,1);
  Q = zeros(m,4);
  T = zeros(m,3);
  for b = 1:m
    fk_helper(b);
  end

end
