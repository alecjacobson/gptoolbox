function [V,F] = execute_triangle(command_line_args, poly_file_name_prefix)
  % execute the util triangle (http://www.cs.cmu.edu/~quake/triangle.html)
  % with the given command line arguements
  %
  % Input:
  % command_line_args: string of command line arguments
  %
  % Outputs:
  %   V  #vertices by 2, list of vertex positions
  %   F  #faces by 3, list of face indices
  %
  % See also: triangle
  %
  warning('THIS FILE IS DEPRECATED. USE TRIANGLE.M INSTEAD');

  % hard code this for my computer for now, need to find out how to find
  % triangle...'which triangle' does not work
  path_to_triangle = '/opt/local/bin/triangle';

  [status, result] =system( ...
    [path_to_triangle ' ' command_line_args ' ' poly_file_name_prefix '.poly']);
  if(status ~= 0 )
    error(result);
  end
  F = readELE([poly_file_name_prefix '.1.ele']);
  V = readNODE([poly_file_name_prefix '.1.node']);
  % Triangle likes to use 1-indexed though .ele reader is 0-indexed
  if(min(F(:)) > 1 && max(F(:)) > size(V,1))
    F = F-1;
  end
  delete( ...
    [poly_file_name_prefix '.1.node'],...
    [poly_file_name_prefix '.1.ele'],...
    [poly_file_name_prefix '.1.poly']);
end
