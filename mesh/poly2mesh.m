function [V,F] = poly2mesh(filename,minimum_angle,maximum_area)
  % POLY2MESH Triangulate interoir of polygon read from .poly  file using
  % Triangle
  %
  %
  % [V,F] = poly2mesh(filename,minimum_angle,maximum_area)
  % 
  % Inputs:
  %  filename  path to .poly file
  %  minimum_angle  minimum angle parameter for Triangle
  %  maximum_area  maximum area parameter for Triangle
  %
  % Outputs:
  %   V  #vertices by 2, list of vertex positions
  %   F  #faces by 3, list of face indices
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: png2objandtga, png2poly, writePOLY, poly2VEH, triangle
  %
  warning('THIS FILE IS DEPRECATED. CALL TRIANGLE.M DIRECTLY INSTEAD');
  basename = regexprep(filename,'\.poly$','');

  minimum_angle_args = '';
  if(exist('minimum_angle'))
    minimum_angle_args = ['q' num2str(minimum_angle)];
  end
  maximum_area_args = '';
  if(exist('maximum_area'))
    maximum_area_args = ['a' num2str(maximum_area)];
  end

  % use triangle to triangulate interior
  command_line_args = ['-p' minimum_angle_args maximum_area_args];
  [V,F] = execute_triangle(command_line_args,basename);
end
