function [w,a] = axisanglebetween(u,v,nan_replacement)
  % AXISANGLEBETWEEN compute the axis and angle representation of the rotation
  % between two vectors
  %
  % [w,a] = axisanglebetween(u,v)
  % [w,a] = axisanglebetween(u,v,nan_replacement) same as above but NaN rows
  % resulting from angle exactly pi are replaced with nan_replacement
  %
  % Inputs:
  %   u  n by dim list of row vectors
  %   v  n by dim list of row vectors
  %   Optional:
  %     nan_replacement vector replacement for NaN rows in w
  % Outputs:
  %   w  n by 3 list of axis vectors, if dim = 2 then this is as if u v had 0s
  %     as z coordinates, NaNs may appear if angle is exactly pi
  %   a  n by 1 list of angles
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %

  assert(all(size(u) == size(v)));
  dim = size(u,2);
  assert((dim == 2) || (dim == 3));
  n = size(u,1);

  u = normalizerow(u);
  v = normalizerow(v);

  if(dim == 2)
    u = [u zeros(n,1)];
    v = [v zeros(n,1)];
  end
  w = normalizerow(cross(u,v,2));
  if(exist('nan_replacement','var'))
    assert(numel(nan_replacement) == 3);
    nan_replacement = reshape(nan_replacement,1,3);
    redo = any(isnan(w),2);
    w(redo,:) = repmat(nan_replacement,[sum(redo) 1]);
  end
  %if(dim == 2)
  %  redo = any(isnan(w),2);
  %  % doesn't matter in 2D just use [0 0 1]
  %  w(redo,:) = repmat([0 0 1],[sum(redo) 1]);
  %end

  a = acos(dot(u,v,2));
  a(normrow(u-v)==0) = 0;
end
