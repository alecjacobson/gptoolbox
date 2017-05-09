function [SI,SO,EI,EO] = offset_curve(P, offset)
  % OFFSET_CURVE  Given a closed curve find an inner and outer offset curve.
  % 
  %[SI,SO,EI,EO] = offset_curve(P, offset)
  %
  % Input:
  %   P  #P by 2 list of curve positions
  %   offset   offset amount
  % Output:
  %   SI  #P by 2 list of inner curve offset positions
  %   SO  #P by 2 list of outer curve offset positions, matches orientation of
  %     input
  %   EI  #P by dim  list of inner curve edges
  %   EO  #P by dim  list of outer curve edges
  %

  % http://www.mathworks.com/matlabcentral/newsreader/view_thread/173494
  x = P(:,1)';
  y = P(:,2)';

  % find slope
  m = diff(x + i*y);
  % use average of slope on each side
  % better result for large spaces between points
  m = [m(1) (m(1:(end-1)) + m(2:end))/2 m(end)];
  
  % calculate offset
  vOff = offset*m*exp(-i*pi/2) ./ abs(m);
  
  % generate output vectors
  IX = fliplr([x-real(vOff)]); 
  IY = fliplr([y-imag(vOff)]);
  OX = [x+real(vOff)];
  OY = [y+imag(vOff)];
  SI = [IX' IY'];
  SO = [OX' OY'];
  EI = [1:numel(IX);2:numel(IX) 1]';
  EO = [1:numel(OX);2:numel(OX) 1]';
end
