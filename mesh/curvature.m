function [kappa,alpha,ev,l] = curvature(P)
  % CURVATURE compute pointwise curvature on a closed piecewise-linear curve P
  %
  % Inputs:
  %   P  #P by 2 list of curve vertices
  % Outputs:
  %   kappa  #P list of curvature values
  %   alpha  #P list of exterior angle values
  %   ev  #P by 2 list of edge vectors, ev(i,:) goes from P(i,:) to P(i+1,:)
  %   l  #P list of edge vector lengths, ev(i,:) goes from P(i,:) to P(i+1,:)
  %

  % edge vectors, evi goes from Pi to Pi+1
  ev = P([2:end 1],:) - P;

  % Exterior angles between consecutive edge vectors
  alpha = atan2( ...
    ev([end 1:end-1],1).*ev(1:end,2) - ev([end 1:end-1],2).*ev(1:end,1), ...
    ev([end 1:end-1],1).*ev(1:end,1) + ev([end 1:end-1],2).*ev(1:end,2) );

  l = sqrt(sum(ev.^2,2));
  kappa = alpha ./ (0.5 * (l([end 1:end-1]) + l));
end
