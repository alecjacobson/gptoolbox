function [q] = project_point_to_line_segment(A,B,p)
  % PROJECT_POINT_TO_LINE_SEGMENT Project a point p onto a line segment AB.
  %
  % Inputs:
  %   A  #A by dim list of source points
  %   B  #A by dim list of dest points
  %   p  point in dim-space
  % Output:
  %   q  closest point to p on the line segment from A to B in dim-space
  %
  warning('Deprecated. Call project_to_lines');

  % vector from A to B
  AB = (B-A);
  % squared distance from A to B
  AB_squared = dot(AB,AB);
  if(AB_squared == 0)
    % A and B are the same point
    q = A;
  else
    % vector from A to p
    Ap = (p-A);
    % from http://stackoverflow.com/questions/849211/
    % Consider the line extending the segment, parameterized as A + t (B - A)
    % We find projection of point p onto the line. 
    % It falls where t = [(p-A) . (B-A)] / |B-A|^2
    t = dot(Ap,AB)/AB_squared;
    if (t < 0.0) 
      % "Before" A on the line, just return A
      q = A;
    else if (t > 1.0) 
      % "After" B on the line, just return B
      q = B;
    else
      % projection lines "inbetween" A and B on the line
      q = A + t * AB;
    end
  end
end
