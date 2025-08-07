function CM = zoe(varargin)
  if nargin < 1
    n = 128;
  else
    n = varargin{1};
  end
  if nargin < 2
    alpha = 0.95;
  else
    alpha = varargin{2};
  end

  t = [0;0.5 - 1e-5;0.5 + 1e-5;1];
  h = ['ffffff'; '489acc'; 'fffae3'; 'ff7424'];
  c = hex2rgb(h);


  CM = interp1(t,c,linspace(0,1,2*n)','linear');
  CM(2:2:end,:) = CM(2:2:end,:) * alpha;
end
