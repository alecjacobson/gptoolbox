function [oP,oC,oI,oT] = spline_offset(P,C,th,tol)
  % SPLINE_OFFSET Offset a cubic Bezier spline by a signed distance, returning
  % the offset curve as another cubic Bezier spline. Each input cubic is
  % replaced by one or more cubics whose end points lie exactly on the offset
  % curve and are tangent to it there; a cubic is recursively split until the
  % distance from its midpoint to the input curve is within tol of |th|.
  %
  % [oP,oC,oI,oT] = spline_offset(P,C,th,tol)
  %
  % Inputs:
  %   P  #P by dim list of control point locations
  %   C  #C by 4 list of indices into P of cubic Bezier curves. Consecutive
  %     curves that share an end point (e.g., C = [1 2 3 4;4 5 6 7]) share the
  %     corresponding point in the output.
  %   th  signed offset distance. Positive offsets to the left of the
  %     direction of travel (the tangent rotated by +90 degrees), negative to
  %     the right.
  %   tol  tolerance on | dist(midpoint of offset piece, input curve) - |th| |
  %     used to decide whether to split a piece further
  % Outputs:
  %   oP  #oP by dim list of control point locations of the offset spline
  %   oC  #oC by 4 list of indices into oP of cubic Bezier curves
  %   oI  #oC list of indices into 1:#C, the input curve each offset cubic came
  %     from
  %   oT  2*#oC list of parameter values; oT(2*j-1:2*j) are the parameters in
  %     C(oI(j),:) at which the j-th offset cubic starts and ends
  %
  % Known issues:
  %   The offset of a curve whose radius of curvature is smaller than |th|
  %   self-intersects. Input curves are assumed to be G1 continuous where they
  %   share end points.
  %
  % Example:
  %   P = [0 0;1 0;2 1;3 1;4 1;5 0;6 0]/6;
  %   C = [1 2 3 4;4 5 6 7];
  %   [oP,oC] = spline_offset(P,C, 0.02,1e-2);
  %   [uP,uC] = spline_offset(P,C,-0.02,1e-2);
  %   plot_spline(P,C);
  %   hold on;
  %   plot_spline(oP,oC);
  %   plot_spline(uP,uC);
  %   hold off;
  %   axis equal;
  %
  % See also: spline_to_poly, plot_spline
  oPT = [];
  oCT = [];
  oI = [];
  oT = [];
  oJ = [];

  for i = 1:size(C,1)
    [oPi,oCi,oiT] = offset_cubic_recursive(P(C(i,:),:),th,tol);
    oCT = [oCT oCi' + size(oPT,2)];
    oPT = [oPT oPi'];
    oJi = nan(size(oPi,1),1);
    oJi([1 end]) = C(i,[1 4]);
    oI = [oI;repmat(i,size(oCi,1),1)];
    oT = [oT;oiT];
    oJ = [oJ;oJi];
  end
  oP = oPT';
  oC = oCT';
  [~,IA,IC] = unique(oJ,'stable');
  oP = oP(IA,:);
  oC = reshape(IC(oC),size(oC));

end

function [oP,oC,oT] = offset_cubic_recursive(C,th,tol)
  oP = offset_cubic(C,th);
  oC = [1 2 3 4];
  oT = [0;1];
  % now check if midpoint is close enough to the offset curve
  om = cubic_eval(oP,0.5);
  dm = sqrt(point_cubic_squared_distance(om,C));
  if abs(dm-abs(th)) < tol
    return;
  end

  [C1,C2] = cubic_split(C,0.5);
  [oP1,oC1,oT1] = offset_cubic_recursive(C1,th,tol);
  oT1 = oT1*0.5;
  [oP2,oC2,oT2] = offset_cubic_recursive(C2,th,tol);
  oT2 = oT2*0.5 + 0.5;
  oP = [oP1(1:end-1,:);oP2];
  oC = [oC1;oC2+size(oP1,1)-1];
  oT = [oT1;oT2];

end

function [oC] = offset_cubic(C,th)
  [T0,T1] = spline_robust_endpoint_tangent_directions(C,[1 2 3 4]);
  N0 = normalizerow([-T0(:,2) T0(:,1)]);
  N1 = normalizerow([-T1(:,2) T1(:,1)]);
  % for a given curve γ(t), the offset curve is o(t) = γ(t) + th*N(t)
  % so we have γ(0) = C(1,:) and γ(1) = C(4,:)
  oC = nan(4,size(C,2));
  oC(1,:) = C(1,:) + th*N0;
  oC(4,:) = C(4,:) + th*N1;
  % we want to find oC(2,:) and oC(3,:) such that the cubic with control points
  % oC is tangent to the offset curve at the endpoints
  %
  % The tangent at t=0 is γ'(0) = ∂γ/∂t(0)
  % So the tangent at t=0 of the offset curve is ∂o/∂t(0) = ∂γ/∂t(0) +
  % th*∂N/∂t(0)
  %
  % if γ is a cubic bezier curve then ∂N/∂t(0) = 
  D0 = cubic_tangent(C,0);
  l0 = normrow(D0);
  D1 = cubic_tangent(C,1);
  l1 = normrow(D1);
  k01 = spline_endpoint_curvature(C,[1 2 3 4]);

  dNdt0 = k01(1)*D0;
  dNdt1 = k01(2)*D1;
  dodt0 = D0 - th*dNdt0;
  dodt1 = D1 - th*dNdt1;

  oC(2,:) = oC(1,:) + dodt0/3;
  oC(3,:) = oC(4,:) - dodt1/3;

  %t = 0
  %fd(@(t) cubic_eval(C,t) + th * normalizerow(cubic_tangent(C,t))*[0 1;-1 0],t)
  %fd(@(t) cubic_eval(oC,t),t)
  
end

