function [C] = slerp(A,B,t)
  % SLERP Spherically interpolate between to unit vectors. If not unit then
  % this will normalize, slerp and then lerp magnitudes.
  %
  % [C] = slerp(A,B,t)
  %
  % Inputs:
  %   A  #A by 3 list of unit vectors
  %   B  #A by 3 list of unit vectors
  %   t  #t list of values between 0 and 1
  % Outputs:
  %   C  #A by 3 list of interpolated unit vectors
  %

  assert(size(A,1)==size(B,1),'#A != #B');
  if numel(t) == 1
    t = repmat(t,size(A,1),1);
  end
  assert(size(A,1)==numel(t),'#A != #T');

  % Original magnitudes
  Amag = sqrt(sum(A.^2,2));
  Bmag = sqrt(sum(B.^2,2));
  % normalize a and b and keep 0's
  A = bsxfun(@rdivide,A,Amag);
  B = bsxfun(@rdivide,B,Bmag);
  Omega = acos( sum(A.*B,2) );
  Omega(abs(Omega-pi)<1e-8) = pi-10*1e-8;
  sinOmega = sin(Omega);
  C = bsxfun(@times,A,sin((1-t).*Omega)./sinOmega) + ...
      bsxfun(@times,B,sin(t.*Omega)./sinOmega);
  % coincident
  co = abs(Omega)<eps;
  % Converge to typical lerp
  if any(co)
    C(co,:) = A(co,:) + bsxfun(@times,t(co),(B(co,:)-A(co,:)));
  end

  % lerp original magnitudes
  C = bsxfun(@times,C,Amag + t.*(Bmag-Amag));

end
