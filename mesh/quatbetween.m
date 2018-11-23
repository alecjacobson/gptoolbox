function q = quatbetween(u,v,varargin)
  % QUATBETWEEN compute the unit quaternion representation of the rotation
  % between two vectors
  %
  % [w,a] = axisanglebetween(u,v)
  % [w,a] = axisanglebetween(u,v,nan_replacement) same as above but NaN rows
  % resulting from angle exactly pi are replaced with nan_replacement
  %
  % Inputs:
  %   u  n by dim list of row vectors
  %   v  n by dim list of row vectors
  % Outputs:
  %   q  n by 4 list of unit quaternions 
  %
  [w,a] = axisanglebetween(u,v,varargin{:});
  q = axisangle2quat(w,a);
end
