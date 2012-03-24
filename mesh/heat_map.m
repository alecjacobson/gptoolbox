function [M,C] = heat_map(varargin)
  % HEAT_MAP create a colormap defined from neg, zero, pos, too_pos fading
  % respectively from blue, white, red and yellow
  %
  % [M,C] = heat_map()
  % [M,C] = heat_map(s)
  % [M,C] = heat_map(s,neg,zero,pos,too_pos)
  %
  % Inputs (optional):
  %   s  number of samples, {128}
  %   neg  negative value, {-1}
  %   zero  zero value, {0}
  %   pos  positive value, {1}
  %   too_pos  too positive value, {2}
  % Outputs:
  %   M  s by 3 list of colormap RGB values between 0 and 1
  %   C  color axes values: [neg too_pos]
  %
  % Example:
  %   % build heat map
  %   [M,C] = heat_map();
  %   % set colormap in current figure
  %   colormap(M);
  %   % set color axis in current figure
  %   caxis(C);
  %
  % See also: colormap
  % 

  % set default input values
  s = 128;
  neg = -1;
  zero = 0;
  pos = 1;
  too_pos = 2;

  % set input if given
  if nargin >= 1
    s = varargin{1};
  end
  if nargin >= 2
    neg = varargin{2};
  end
  if nargin >= 3
    zero = varargin{3};
  end
  if nargin >= 4
    pos = varargin{4};
  end
  if nargin >= 5
    too_pos = varargin{5};
  end
  if nargin >= 6
    error('Too many inputs');
  end

  assert(neg <= zero);
  assert(zero <= pos);
  assert(pos <= too_pos);
  assert(neg < too_pos);

  % range of remaining values
  r = too_pos - neg;
  % steps per interval
  I = zeros(3,1);
  % remaining steps
  rem = s;
  I(1) = round(((zero-neg)/r)*rem);
  rem = rem - I(1);
  r = too_pos - zero;
  I(2) = round(((pos-zero)/r)*rem);
  rem = rem - I(2);
  I(3) = rem;

  assert(sum(I) == s);

  NEG = [0 0 1];
  ZERO = [1 1 1];
  POS = [1 0 0];
  TOO_POS = [1 1 0];

  neg_zero = ...
    repmat(linspace(1,0,I(1)+1)',1,3).*repmat(  NEG,I(1)+1,1) + ...
    repmat(linspace(0,1,I(1)+1)',1,3).*repmat( ZERO,I(1)+1,1);
  zero_pos = ...
    repmat(linspace(1,0,I(2)+1)',1,3).*repmat( ZERO,I(2)+1,1) + ...
    repmat(linspace(0,1,I(2)+1)',1,3).*repmat(  POS,I(2)+1,1);
  pos_too_pos = ...
    repmat(linspace(1,0,I(3))',1,3).*repmat(    POS,I(3),1) + ...
    repmat(linspace(0,1,I(3))',1,3).*repmat(TOO_POS,I(3),1);

  M = [ ...
    neg_zero(1:(end-1),:); ...
    zero_pos(1:(end-1),:); ...
    pos_too_pos(1:end,:)]; 

  assert(size(M,1) == s);

  C = [neg too_pos];

end
