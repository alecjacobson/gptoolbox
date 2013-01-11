function writePOLY(varargin)
  % WRITEPOLY **obselete** 
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: writePOLY_tetgen, writePOLY_triangle, writePOLY_pyramid
  %
  error(['This is obselete. ' ...
    'Call writePOLY_tetgen.m, *_triangle.m, or *_pyramid.m directly']);
  % dummy calls so that writePOLY_*s get included when zipping dependencies
  writePOLY_tetgen();
  writePOLY_triangle();
  writePOLY_pyramid();
end
