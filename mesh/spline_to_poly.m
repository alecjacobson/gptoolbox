function [V,E,I,K,T] = spline_to_poly(P,C,tol)
  % SPLINE_TO_POLY Evaluate a cubic Bezier spline as a polyline where each
  % segment corresponds to a locally flat segment of the curve up to given
  % tolerance
  % 
  % [V,E] = spline_to_poly(P,C,tol)
  % 
  % Inputs:
  %   P  #P by dim list of control point locations
  %   C  #C by 4 list of indices into P of cubic Bezier curves
  %   tol  tolerance 
  % Outputs:
  %   V  #V by dim list of vertex locations
  %   E  #E by dim list of edge indices into V
  %   I  #E list of indices into 1:#C
  %   K  #V list of indices into 1:#C
  %   T  #V list of parametric locations. T(i) is the parametric location of
  %     V(i,:) in the cubic Bezier curve C(K(i),:)
  %
  % Example:
  %  rng(4);
  %  P = randn(7,2)*400; P(5,:) = P(4,:) + abs(randn(1))*(P(4,:)-P(3,:));
  %  C = [1 2 3 4; 4 5 6 7];
  %  [V,E,I,K,T] = spline_to_poly(P,C,10);
  %  EI = (1:size(E,1))';
  %  B = rand(size(E,1),1)*0.25+0.5;
  %  J = I(EI);
  %  T2 = T(E(EI,2));
  %  T2(J ~= K(E(EI,2))) = 1;
  %  S = T(E(EI,1)) + B(:,2).*(T2-T(E(EI,1)));
  %  X = cell2mat(arrayfun(@(i) cubic_eval(P(C(J(i),:),:),S(i)),(1:size(B,1))','UniformOutput',false));

  % use transpose for amortized concatenation...
  V = P';
  T = nan(size(P,1),1);
  K = zeros(size(P,1),1);
  T(C(:,4)) = 1;
  K(C(:,4)) = 1:size(C,1);
  T(C(:,1)) = 0;
  K(C(:,1)) = 1:size(C,1);

  E = [];
  I = [];
  % consider each cubic
  for c = 1:size(C,1)
    [Pc,Tc] = cubic_flat_eval(P(C(c,:),:),tol);
    J = [C(c,1) size(V,2)+(1:size(Pc,1)-2) C(c,4)];
    Ec = [J(1:end-1);J(2:end)]';
    % use transpose for amortized concatenation...
    E = [E,                    Ec'];
    V = [V,         Pc(2:end-1,:)'];
    T = [T;Tc(2:end-1,:)];
    K = [K;repmat(c,size(Pc,1)-2,1)];
    I = [I;repmat(c,size(Ec,1),1)];
  end
  V = V';
  E = E';
  [V,~,keep,E] = remove_unreferenced(V,E);
  T = T(keep);
  K = K(keep);
  assert(~any(isnan(T)));

end
