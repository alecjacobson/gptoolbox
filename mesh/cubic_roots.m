function Z = cubic_roots(C)
  % Z = cubic_roots(C)
  %
  % https://en.wikipedia.org/wiki/Cubic_function#General_formula
  % ξ^[2,1,0]
  a = C(:,1);
  b = C(:,2);
  c = C(:,3);
  d = C(:,4);

  % Fails if H==0
  xi = [1 -0.5+0.5*sqrt(3)*i -0.5-0.5*sqrt(3)*i ];
  D = 18.*a.*b.*c.*d - 4.*b.^3.*d + b.^2.*c.^2 -4.*a.*c.^3 -27.*a.^2.*d.^2;
  D0 = b.^2 - 3*a.*c;
  D1 = 2.*b.^3 - 9.*a.*b.*c + 27.*a.^2.*d;
  Z = nan(size(C,1),1)+i;
  K = sqrt(D1.^2 - 4.*D0.^3);
  H = ((D1 + K)/2).^(1/3);
  % This seems like a "hot fix" but I don't have documentation why I don't need
  % to worry about the cases below
  HP = ((D1 - K)/2).^(1/3);
  H(abs(H)<eps) = HP(abs(H)<eps);
  Z = -1./(3.*a).*(b + H.*xi + D0./(xi.*H));
  %ZDD0 = repmat(-b./(3.*a),1,3)
  %ZD = [ ...
  %  repmat((9.*a.*d - b.*c)./(2.*D0),1,2) ...
  %  (4.*a.*b.*c - 9.*a.^2.*d - b.^3)./(a.*D0)]
  %Z( abs(D)<eps,:) = ZD( abs(D)<eps,:);
  %Z( abs(D)<eps & abs(D0)<eps,:) = ZDD0( abs(D)<eps & abs(D0)<eps,:);

  %% % [x] = cubic (a,b,c,d)
  %%
  %%   Gives the roots of the cubic equation
  %%         ax^3 + bx^2 + cx + d = 0    (a <> 0 !!)
  %%   by Nickalls's method: R. W. D. Nickalls, ``A New Approach to
  %%   solving the cubic: Cardan's solution revealed,''
  %%   The Mathematical Gazette, 77(480)354-359, 1993.
  %%   dicknickalls@compuserve.com
  %%  Herman Bruyninckx 10 DEC 1996, 19 MAY 1999
  %%  Herman.Bruyninckx@mech.kuleuven.ac.be 
  %%  Dept. Mechanical Eng., Div. PMA, Katholieke Universiteit Leuven, Belgium
  %%  <http://www.mech.kuleuven.ac.be/~bruyninc>
  %%
  %% THIS SOFTWARE COMES WITHOUT GUARANTEE.
  %% https://people.mech.kuleuven.be/~bruyninc/matlab/cubic.m
  
  %x = nan(size(C,1),3);
  %
  %xN = -b./3./a;
  %yN = d + xN .* (c + xN .* (b + a.*xN));
  %
  %two_a    = 2*a;
  %delta_sq = (b.*b-3.*a.*c)./(9.*a.*a);
  %h_sq     = two_a .* two_a .* delta_sq.^3;
  %dis      = yN.*yN - h_sq;
  %onethird      = 1/3;



  %% three real roots (two or three equal):
  %delta = (yN./two_a).^onethird;
  %Z = [xN + delta  xN + delta  xN - 2*delta];

  %% three distinct real roots:
  %theta = acos(-yN./sqrt(h_sq))./3;
  %delta = sqrt(delta_sq);
  %two_d = 2*delta;
  %twop3 = 2*pi/3;
  %x_dis_lt_m_eps = [ ...
  %  xN + two_d.*cos(theta) ...
  %  xN + two_d.*cos(twop3-theta) ...
  %  xN + two_d.*cos(twop3+theta)];
  %dis_lt_m_eps = dis < -eps;
  %Z(dis_lt_m_eps,:) = x_dis_lt_m_eps(dis_lt_m_eps,:);
  %
  %% one real root:
  %dis_sqrt = sqrt(dis);
  %r_p  = yN - dis_sqrt;
  %r_q  = yN + dis_sqrt;
  %p    = -sign(r_p) .* ( sign(r_p).*r_p./two_a ).^onethird;
  %q    = -sign(r_q) .* ( sign(r_q).*r_q./two_a ).^onethird;
  %x2 = xN + p * exp(2*pi*i/3) + q * exp(-2*pi*i/3);
  %x_dis_ge_eps = [ xN+p+q x2 conj(x2)]
  %dis_ge_eps = dis >= eps;
  %Z(dis_ge_eps,:) = x_dis_ge_eps(dis_ge_eps,:);


end
