function [d,pair] = hausdorff(VA,FA,VB,FB)
  % HAUSDORFF compute the Hausdorff distance between mesh (VA,FA) and mesh
  % (VB,FB). This is the 
  %
  % d(A,B) = max ( max min d(a,b) , max min d(b,a) )
  %                a∈A b∈B          b∈B a∈A
  %
  % Inputs:
  %   VA  #VA by 3 list of vertex positions
  %   FA  #FA by 3 list of face indices into VA
  %   VB  #VB by 3 list of vertex positions
  %   FB  #FB by 3 list of face indices into VB
  % Outputs:
  %   d  hausdorff distance
  %   pair  2 by 3 list of "determiner points" so that pair(1,:) is from A and
  %     pair(2,:) is from B
  [sqr_DBA,~,CBA] = point_mesh_squared_distance(VB,VA,FA);  
  [max_sqrDBA,ind] = max(sqr_DBA);
  max_CBA = [CBA(ind,:);VB(ind,:)];
  [sqr_DAB,~,CAB] = point_mesh_squared_distance(VA,VB,FB);  
  [max_sqrDAB,ind] = max(sqr_DAB);
  max_CAB = [VA(ind,:);CAB(ind,:)];
  if max_sqrDBA > max_sqrDAB
    d = sqrt(max_sqrDBA);
    pair = max_CBA;
  else
    d = sqrt(max_sqrDAB);
    pair = max_CAB;
  end
end
