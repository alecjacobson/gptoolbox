function [V,E,cid] = get_pencil_curves(tol)
  % GET_PENCIL_CURVES Get many pencil curves
  %
  % [V,E,cid] = get_pencil_curves(tol)
  %
  % Inputs:
  %   tol  dpsimplify tolerance
  % Outputs:
  %   V  #V by 2 list of vertices
  %   E  #E by 2 list of edge indices into V
  %   cid list of indices into V indicating last vertex in each curve
  %
  % See also: get_pencil_curve
  %
  V = [];
  E = [];
  cid = [];
  while true
    hold on;
    V1 = get_pencil_curve();
    hold off;
    V1 = dpsimplify(V1,tol);
    E1 = [1:size(V1,1)-1;2:size(V1,1)]';
    E = [E;size(V,1)+E1];
    V = [V;V1];
    cid = [cid;size(V,1)];
    if input('Continue? y/n','s') == 'n'
      return
    end
  end
end
