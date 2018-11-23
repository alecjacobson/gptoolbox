function G = grad_intrinsic(l,F)
  % GRAD_INTRINSIC Construct an intrinsic gradient operator.
  %
  % Inputs:
  %  l  #F by 3 list of edge lengths
  %  F  #F by 3 list of triangle indices into some vertex list V
  % Outputs:
  %  G  #F*2 by #V gradient matrix: G=[Gx;Gy] where x runs along the 23 edge and
  %    y runs in the counter-clockwise 90° rotation.
  %
  
  x = (l(:,2).^2-l(:,1).^2-l(:,3).^2)./(-2.*l(:,1));
  y = sqrt(l(:,3).^2 - x.^2);
  n = max(F(:));
  m = size(F,1);
  Z = zeros(m,1);
  V2 = [x y;Z Z;l(:,1) Z];
  F2 = reshape(1:3*m,m,3);
  G2 = grad(V2,F2);
  P = sparse(F2,F,1,m*3,n);
  G = G2*P;

end
