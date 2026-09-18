  % [sqrD,S,K] = point_cubic_squared_distance(Q,C)
  %
  % Inputs:
  %   Q   #Q by dim list of query points
  %   C   4 by dim list of control points
  % Outputs:
  %   sqrD  #Q list of squared distances
  %   S  #Q list of parameter value on cubic spline
  %   K  #Q by dim list of closest points on corresponding cubic spline
