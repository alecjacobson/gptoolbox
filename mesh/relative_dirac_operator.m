function DR = relative_dirac_operator(V,F)
  % RELATIVE DIRAC OPERATOR Construct the discrete quaternionic relative dirac
  % operator for a 3d triangle mesh (as described by "A Dirac Operator for
  % Extrinsic Shape Analysis" [Liu, Jacobson, and Crane. 2017]).
  % 
  % Inputs:
  %   V   #V by 3 list of vertex positions
  %   F   #F by 3 list of triangle indices into V
  % Outputs:
  %   DR  4*#F by 4*#V sparse rectangular relative Dirac matrix. Each 4x4
  %     sublock represents a quaternion (w,x,y,z)
  %

  A = doublearea(V,F);
  N = per_vertex_normals(V,F);
  dN = [zeros(numel(F),1) N(F(:,[2 3 1]),:) - N(F(:,[3 1 2]),:)];
  % [ w -x -y -z, x  w -z  y, y  z  w -x, z -y  x  w];
  Q = [ ...
    1 0 0 0;0 -1  0 0;0 0 -1  0;0  0 0 -1; ...
    0 1 0 0;1  0  0 0;0 0  0 -1;0  0 1  0; ...
    0 0 1 0;0  0  0 1;1 0  0  0;0 -1 0  0; ...
    0 0 0 1;0  0 -1 0;0 1  0  0;1  0 0  0]';
  % face-quat indices 
  II = repmat(repmat((0:size(F,1)-1)'*4,1,4) + repmat((1:4),size(F,1),1),3,4);
  JJ = (repmat(F(:),1,4*4)-1)*4 + repmat(reshape(repmat(1:4,4,1),1,[]),size(F(:),1),1);
  DR = sparse(II,JJ,-dN*Q./repmat([A;A;A],1,16),4*size(F,1),4*size(V,1));
end
