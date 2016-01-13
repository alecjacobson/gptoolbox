function [E] = sharp_edges(V,F,varargin)
  % SHARP_EDGES Given a mesh, compute sharp edges.
  %
  % [E] = sharp_edges(V,F)
  % [E] = sharp_edges(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions 
  %   F  #F by 3 list of mesh triangle indices into V
  %   Optional:
  %     'Angle'  followed by dihedral angle considered sharp.
  % Outputs:
  %   E  #E by 2 list of edge indices into V {pi*0.11}
  %

  angle = pi*0.11;

  % sharp edges
  [A,C] = adjacency_dihedral_angle_matrix(V,F);
  %% This is much much slower
  %A(1&A) = abs(A(1&A)-pi)>pi*0.11;
  [AI,AJ,AV] = find(A);
  keep = abs(AV-pi)>(angle) & ~isnan(AV);
  A = sparse(AI(keep),AJ(keep),1,size(A,1),size(A,2));
  [CI,~,CV] = find(C.*A);
  II = [CI+mod(CV,3)*size(F,1) CI+mod(CV+1,3)*size(F,1)];
  E = F(II);

end
