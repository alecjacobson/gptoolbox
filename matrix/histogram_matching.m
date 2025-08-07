function Z = histogram_matching(X,Y,sigma)
  % Z = histogram_matching(X,Y,sigma)
  %
  % Inputs:
  %   X  #X histogram values (assumed to be on [0,1])
  %   Y  #Y values to be matched to (assumed to be on [0,1])
  %   sigma  optional smoothing parameter (default: 0)
  % Outputs:
  %   Z  #X values matched to Y
  % 
  % Villani p. 73

  % Z = H_Y⁻¹(H_X(X))
  dims = size(X);
  X = X(:);
  Y = Y(:);

  [x,I] = sort(X);
  % H_X(x) = t_x
  t_x = linspace(0, 1, numel(x))';
  y = unique(Y);
  % H_Y(y) = t_y
  t_y = linspace(0, 1, numel(y))';

  y0 = y;


  % z = H_Y⁻¹(H_X(x))
  z = interp1(t_y, y, t_x, 'linear', 'extrap');

  dx = z-x;

  % This seems reasonable
  if nargin > 2 && sigma > 0
    E = [1:numel(dx)-1; 2:numel(dx)]';
    L = cotmatrix(t_x,E);
    M = massmatrix(t_x,E);
    b = [1;numel(dx)];
    bc = [dx(1);dx(end)];
    w = sigma^2;
    dx = min_quad_with_fixed(-0.5*w*L+0.5*M,-M*dx,b,bc);
  end

  Z(I) = x+dx;
  
  %clf;
  %hold on;
  %plot(x, t_x, 'r');
  %plot(y0, t_y, ':b');
  %plot(y, t_y, 'b');
  %plot(z, t_x, 'g');
  %plot(z-x,t_x,'k');
  %v = unique(Z);
  %% H_Y(y) = t_y
  %t_v = linspace(0, 1, numel(v))';
  %plot(v, t_v, 'm');
  %hold off;
  %error

  Z = reshape(Z, dims);
end
