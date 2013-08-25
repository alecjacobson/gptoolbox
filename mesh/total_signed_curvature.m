function [tsc,alpha,V] = total_signed_curvature(P)
  % TOTAL_SIGNED_CURVATURE  Total discrete signed curvature of a given closed
  % planar curve P (see http://ddg.cs.columbia.edu/SIGGRAPH05/Didactic.pdf)
  %
  % [tsc,alpha] = total_signed_curvature(P)
  %
  % Inputs:
  %   P  #P by 2 a list of positions of a closed piecewise-linear curve
  % Outputs:
  %   tsc  total signed curvature scalar, should be 2*pi*k where k∈0,1,2,...
  %   alpha  #P by 1 list of exterior angles
  %
  %

  % edge vectors
  V = P([2:end 1],:)-P;

  % Exterior angles between consecutive edge vectors
  alpha = atan2( ...
    V(1:end,1).*V([2:end 1],2) - V(1:end,2).*V([2:end 1],1), ...
    V(1:end,1).*V([2:end 1],1) + V(1:end,2).*V([2:end 1],2) );
  
  % total signed curvature
  tsc = sum(alpha);

end
