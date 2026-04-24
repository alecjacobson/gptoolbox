function max_k = cubic_max_curvature_single_interval(C)

%% Power basis
d = C(1,:);
c = 3*(C(2,:) - C(1,:));
b = 3*(C(3,:) - 2*C(2,:) + C(1,:));
a = C(4,:) - 3*C(3,:) + 3*C(2,:) - C(1,:);

ax = a(1); ay = a(2);
bx = b(1); by = b(2);
cx = c(1); cy = c(2);

px = [3*ax 2*bx cx];
py = [3*ay 2*by cy];

qx = [6*ax 2*bx];
qy = [6*ay 2*by];

det_poly = conv(px,qy) - conv(py,qx);
speed2 = conv(px,px) + conv(py,py);

detp = polyder(det_poly);
speed2p = polyder(speed2);

A = conv(detp,speed2);
B = 1.5*conv(det_poly,speed2p);

n = max(length(A),length(B));
A = [zeros(1,n-length(A)) A];
B = [zeros(1,n-length(B)) B];

crit_poly = A - B;

R = fast_roots(crit_poly,0,1);
R = R(~isnan(R));

R = R(:);
candidates = [0;1;R];

dx  = polyval(px,candidates);
dy  = polyval(py,candidates);
ddx = polyval(qx,candidates);
ddy = polyval(qy,candidates);

num = dx.*ddy - dy.*ddx;
den = (dx.^2 + dy.^2).^(3/2);

k = abs(num ./ den);
k(~isfinite(k)) = 0;

max_k = max(k);
end
