function P = cubic_eval(C,t)
  % CUBIC_EVAL Evaluate a cubic Bezier curve.
  %
  % P = cubic_eval(C,t);
  %
  % Inputs:
  %   C  4 by dim list of control points
  %   t  #t list of evaluation parameters
  % Outputs:
  %   P  #t by dim list of evaluated points
  %
  % See also: readSVG_cubics, cubic_split
  %
  t = reshape(t,numel(t),1);
  %B = bernstein_eval([0 1 2 3],3,t);
  %P = B*C;
  P =  ...
    1*(1-t).^3.*t.^0.*C(1,:) +  ...
    3*(1-t).^2.*t.^1.*C(2,:) +  ...
    3*(1-t).^1.*t.^2.*C(3,:) +  ...
    1*(1-t).^0.*t.^3.*C(4,:);
end
