function [dS,Q] = retarget(C,E,P,dQ,SC)
  % RETARGET Retarget a set of relative bone transformations (dQ) from one
  % skeleton (C,E,P) to another skeleton (SC,E,P)
  %
  % Inputs:
  %   C  #C by 3 list of joint positions
  %   E  #E by 2 list of bone edge indices into C
  %   P  #E list of parent bone indices into E (0 means root)
  %   dQ  #E by 4 list of relative bone transformations as quaternions
  %   SC  #C by 3 list of joint positions
  % Outputs:
  %   dS  #E by 4 list of relative bone transformations as quaternions
  %   Q  #E by 4 list of "change of basis" transformations mapping each edge in
  %     (C,E,P) to its corresponding edge in (SC,E,P). Such that: 
  %     dS = quatmultiply(quatconj(Q),quatmultiply(dQ,Q));
  %  

  % Find rotation that takes other bone frame to this bone frame
  EV = normalizerow(C(E(:,2),:)-C(E(:,1),:));
  SEV = normalizerow(SC(E(:,2),:)-SC(E(:,1),:));
  [W,A] = axisanglebetween(SEV,EV,[0 1 0]);
  Q = axisangle2quat(W,A);

  % Using the FK hierarchy **does** produce a different change of basis. I
  % guess this corresponds to the parallel transport along the FK chain, which
  % seems like a reasonable thing to do.
  Q = zeros(size(E,1),4);
  seen = false(size(E,1),1);
  function Qq = change(q)
    if ~seen(q)
      seen(q) = true;
      Qp = [1 0 0 0];
      if P(q) > 0
        change(P(q));
        Qp = Q(P(q),:);
      end
      % q e q' = s
      % q = d p 
      % (dp)e(dp)' = s
      %  d pep' d' = s
      s = SEV(q,:);
      e = EV(q,:);
      [W,A] = axisanglebetween(s,quatrotate(Qp,e),[0 1 0]);
      Qd = axisangle2quat(W,A);
      Q(q,:) = quatmultiply(Qd,Qp);
      [W,A] = axisanglebetween(s,e);
    end
  end
  for q = 1:size(Q,1)
    change(q);
  end

  dS = quatmultiply(quatconj(Q),quatmultiply(dQ,Q));
end
