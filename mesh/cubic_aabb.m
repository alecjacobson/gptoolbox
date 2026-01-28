function [C1,C2] = cubic_aabb(C)
  % [C1,C2] = cubic_aabb(C)
  %
  % Inputs:
  %   C  4 by dim list of control points
  % Outputs:
  %   C1  1 by dim of minimum values
  %   C2  1 by dim of maximum values
  %
  function [ts,A,B,C] = extrema_ts(U)
    A = 9*U(2) - 3*U(1) - 9*U(3) + 3*U(4);
    B = 6*U(1) - 12*U(2) + 6*U(3);
    C = -3*U(1) + 3*U(2);
    ts = roots([A B C]);
    ts = ts(isreal(ts) & ts>0 & ts<1);
  end

  C1 = inf(1,size(C,2));
  C2 = -inf(1,size(C,2));
  for c = 1:size(C,2)
    ts = [extrema_ts(C(:,c))];
    Qc = cubic_eval(C(:,c),ts);
    C1(1,c) = min([C1(:,c);Qc;C(1,c);C(4,c)]);
    C2(1,c) = max([C2(:,c);Qc;C(1,c);C(4,c)]);
  end
end

