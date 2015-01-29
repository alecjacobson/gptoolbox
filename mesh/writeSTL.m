function writeSTL(filename, V,F)
  % WRITESTL writes an STL file with vertex/face information. Uses single
  % precision.
  %
  % writeSTL(filename,V,F)
  %
  % Input:
  %  filename  path to .stl file
  %  V  #V by 3 list of vertices (Default units are in mm)
  %  F  #F by 3 list of triangle indices
  %
  stlwrite(filename,F,V);
end

