function P = quadratic_eval(Q,t)
  % QUADRATIC_EVAL Evaluate a cubic Bezier curve.
  %
  % P = quadratic_eval(C,t);
  %
  % Inputs:
  %   Q  3 by dim list of control points
  %   t  #t list of evaluation parameters
  % Outputs:
  %   P  #t by dim list of evaluated points
  %
  % See also: readSVG_cubics, cubic_split, quadratic_eval
  %
  t = reshape(t,numel(t),1);
  P = ...
    (1-t).^2   .*Q(1,:) +  ...
    2.*(1-t).*t   .*Q(2,:) +  ...
           t.^2.*Q(3,:);
end
