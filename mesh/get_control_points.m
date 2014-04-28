function [C,flag] = get_control_points(V,F)
  % GET_CONTROL_POINTS displays a mesh and prompts user to select positions of
  % controls points with the mouse. Press enter when finished or escape to
  % cancel
  %
  % [C] = get_control_points(V,F)
  %
  % Inputs:
  %   V  #V by dim list of control points
  %   F  #F by 3 list of face indices
  % Outputs:
  %   C  #C by 2 list of control point positions, if escaped is pressed then
  %     this will be empty
  %   flag  false only if getpts was cancelled
  % 
  % See also: getpts
  %

  V(:,end+1:3) = 0;
  trisurf(F,V(:,1),V(:,2),V(:,3));
  view(2);
  axis equal;

  Cx = zeros(0,1);
  Cy = zeros(0,1);
  flag = false;
  try
    [Cx,Cy] = getpts();
  flag = true;
  catch e
      e.identifier
    switch e.identifier
    case {'MATLAB:license:checkouterror','MATLAB:UndefinedFunction'}
      % probably don't have image processing toolbox
      [Cx,Cy] = ginput;
    end
  end

  C = [Cx Cy];
end
