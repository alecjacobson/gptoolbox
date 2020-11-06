function [P,C,I] = cubic_subdivide(B,T)
  % CUBIC_SUBDIVIDE Given a cubic Bezier curve C split at all given paramter
  % values in T, producing a spline with shared control points and T+1
  % sub-curves.
  %
  % [P,C] = cubic_subdivide(B,T)
  %
  % Inputs:
  %   B  #4 by dim list of control point locations
  %   T  #T list of parameters between [0,1] 
  % Outputs:
  %   P  #P by dim list of control points, P(1:2,:) == B([1 4],:)
  %   C  3(#T+1)+1 by 4 list of control point indices into P
  %   I  #P index into [B; [all new points] ]
  %
  % See also: cubic_split, cubic_eval
  P = B;
  C = [1 2 3 4];
  T = sort(T);
  ca = 1;
  Ca = C(ca,:);
  for ti = 1:numel(T)
    [C1,C2] = cubic_split(P(Ca,:),T(ti));
    Ca1 = [Ca(1) size(P,1)+(1:3)];
    Ca = [size(P,1)+(3:5) Ca(4)];
    P = [P;C1(2:4,:);C2(2:3,:)];
    C(ca,:) = Ca1;
    C = [C;Ca];
    ca = size(C,1);
    % adjust parameter
    T(ti+1:end) = (T(ti+1:end) - T(ti))./(1-T(ti));
  end
  [P,~,I,C] = remove_unreferenced(P,C);
  C = reshape(C,[],4);
end
