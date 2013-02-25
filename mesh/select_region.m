function [S,P] = select_region(t,P)
  % SELECT_REGION select a region of a mesh viewed in a matlab trisurf plot.
  % Assume that the user has already selected a viewing angle and has a handle
  % to the trisurf 'patch'
  %
  % S = select_region(t)
  % [S,P] = select_region(t,P)
  % 
  % Inputs:
  %  t  handle to trisurf object
  %  Optional:
  %    P  polygon from previous call
  % Ouputs:
  %  C  #V list of logicals, true if in selection, false if not
  %  P  polygon to recall this function without prompting user
  %  
  % Example:
  %  % Display mesh (V,F)
  %  t = trisurf(F,V(:,1),V(:,2),V(:,3));axis equal
  %  % User should now have picked a nice view
  %  S = select_region(t);
  %  % Update display to be colored based on selection
  %  set(t,'CData',S*1);
  %

  % get view of axis of handle t
  [az,el] = view(get(t,'Parent'));
  % get figure of t
  of = get(get(t,'Parent'),'Parent');

  % get vertices
  V = get(t,'Vertices');
  % get vertices
  F = get(t,'Faces');

  % rotate by 90 degrees
  Rx90 = [1 0 0; 0 cos(-pi/2) -sin(-pi/2);0 sin(-pi/2) cos(-pi/2)];
  % rotatation according to elevation from view of t
  Rx = [1 0 0; 0 cos(el*pi/180) -sin(el*pi/180);0 sin(el*pi/180) cos(el*pi/180)];
  % rotation according to azimuth of view of t
  Rz = [cos(-az*pi/180) -sin(-az*pi/180) 0; sin(-az*pi/180) cos(-az*pi/180) 0;0 0 1];

  % rotation matrix that will rotate mesh to 
  R = (Rx90*Rx*Rz);
  R = R(1:2,:);
  view_V = [V] * R';
  % use zero for z values
  view_V(:,3) = 0;

  
  if ~exist('P','var') || isempty(P)
  % make a new figure
  f = figure;
  set(f,'Position',get(of,'Position'));
  
  trisurf(F,view_V(:,1),view_V(:,2),view_V(:,3));
  axis equal;
  view(2)
  % get window size
  A = axis;
  A = reshape(axis,2,[]);
  % grow by 100% about center
  dA = [A(1,:)-A(2,:);A(2,:)-A(1,:)]/2;
  axis(reshape(A+dA,1,prod(size(A))));
  
  title('Make your selection here in this figure window:');
  % ask user for a polygon
  h = impoly;
  % Get positions of polygon vertices, before closing figure
  P = getPosition(h);
  % close selection figure
  close(f);
  end

  % In projection (screen space), determine mesh vertices are in selection or
  % not
  S = inpolygon(view_V(:,1),view_V(:,2),P(:,1),P(:,2));
  % convert to double
  S = 1.0*S;

end
