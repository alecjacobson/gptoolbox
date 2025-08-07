function T = cubic_tangent(C,t)
  % T = cubic_tangent(C,t)
  % 
  % Inputs:
  %   C  4 by dim list of control points
  %   t  #t list of parameter values
  % Outputs:
  %   T  #t by dim list of tangent vectors
  t = reshape(t,numel(t),1);
  T = ...
    3*(1-t).^2.*t.^0.*(C(2,:)-C(1,:)) + ...
    6*(1-t).^1.*t.^1.*(C(3,:)-C(2,:)) + ...
    3*(1-t).^0.*t.^2.*(C(4,:)-C(3,:));
end
