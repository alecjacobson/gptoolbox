function [P,T] = cubic_flat_eval(C,tol)
  % CUBIC_FLAT_EVAL Recursively subdivide a cubic Bezier curve until each
  % segment spans a region of the curve that is locally flat up to a given
  % tolerance (i.e., this computes an adaptive refinement).
  %
  % [P,T] = cubic_flat_eval(C,tol)
  %
  % Inputs:
  %   C  4 by dim list of control points
  % Outputs:
  %   P  #P by dim list of evaluated points
  %   T  #T list of corresponding parameter values
  % 
  % See also: cubic_eval, cubic_is_flat, cubic_split
  %
  if cubic_is_flat(C,tol)
    P = C([1 4],:);
    T = [0;1];
  else
    [C1,C2] = cubic_split(C,0.5);
    [P1,T1] = cubic_flat_eval(C1,tol);
    [P2,T2] = cubic_flat_eval(C2,tol);
    P = [P1;P2(2:end,:)];
    T = [T1*0.5;0.5+0.5*T2(2:end)];
  end
end
