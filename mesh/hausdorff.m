function [d,pair] = hausdorff(VA,FA,VB,FB)
  % HAUSDORFF compute the Hausdorff distance between mesh (VA,FA) and mesh
  % (VB,FB). This is the 
  %
  % d(A,B) = max ( max min d(a,b) , max min d(b,a) )
  %                a∈A b∈B          b∈B a∈A
  %
  %  Known issue: This is only computing max(min(va,B),min(vb,A)). This is
  %  better than max(min(va,Vb),min(vb,Va)). This (at least) is missing
  %  "edge-edge" cases like the distance between the two different
  %  triangulations of a non-planar quad in 3D. Even simpler, consider the
  %  Hausdorff distance between the non-convex, block letter V polygon (with 7
  %  vertices) in 2D and its convex hull. The Hausdorff distance is defined by
  %  the midpoint in the middle of the segment across the concavity and some
  %  non-vertex point _on the edge_ of the V.
  % 
  %  Inputs:
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
