function [T,TI,E] = polygons_to_triangles(PI,PC)
  % POLYGONS_TO_TRIANGLES Given a polygonal mesh, replace each polygon with a
  % fan of triangles.
  %
  % [T,TI,E] = polygons_to_triangles(PI,PC)
  % 
  % Inputs:
  %   PI  #PI stream of polygon indices into rows of some V
  %   PC  #polys+1 list of cumulative sum of dual face valences
  % Outputs:
  %   T  #T by 3 list of triangle indices into rows of V
  %   TI  #TI list of indices into input faces (1:#polys)
  %   E  #E by 2 list of edge indices into rows of V
  %
  % See also: dual
  I1 = repelem(PI(PC(1:end-1)+1),diff(PC)-2);
  I2 = PI;
  % remove first and last entries
  I2([PC(1:end-1)+1;PC(2:end)]) = [];
  I3 = PI;
  % remove first two entries
  I3([PC(1:end-1)+1;PC(1:end-1)+2]) = [];
  T = [I1 I2 I3];
  TI = repelem(1:numel(PC)-1,diff(PC)-2);
  if nargout>2
    I1 = PI;
    % order index per polygon
    Pj = [(1:numel(PI))' - repelem(PC(1:end-1),diff(PC))];
    % increment by one mod size
    Pj = mod(Pj,repelem(diff(PC),diff(PC)))+1;
    I2 = PI(repelem(PC(1:end-1),diff(PC)) + Pj);
    E = [I1 I2];
  end
end
