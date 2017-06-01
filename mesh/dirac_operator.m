function D = dirac_operator(V,F)
  % DIRAC_OPERATOR Construct the discrete quaternionic dirac operator for a 3d
  % triangle mesh (as described by "Spin Transformations of Discrete Surfaces"
  % [Crane et al. 2011]).
  % 
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  % Outputs:
  %   D  4*#F by 4*#V sparse rectangular Dirac matrix. Each 4x4 sublock
  %     represents a quaternion (w,x,y,z)
  %

  A = doublearea(V,F);
  EV = [zeros(numel(F),1) V(F(:,[2 3 1]),:) - V(F(:,[3 1 2]),:)];
  % [ w -x -y -z, x  w -z  y, y  z  w -x, z -y  x  w];
  Q = [ ...
    1 0 0 0;0 -1  0 0;0 0 -1  0;0  0 0 -1; ...
    0 1 0 0;1  0  0 0;0 0  0 -1;0  0 1  0; ...
    0 0 1 0;0  0  0 1;1 0  0  0;0 -1 0  0; ...
    0 0 0 1;0  0 -1 0;0 1  0  0;1  0 0  0]';
  % face-quat indices 
  II = repmat(repmat((0:size(F,1)-1)'*4 +(1:4),1,4),3,1);
  JJ = (repmat(F(:),1,4*4)-1)*4+reshape(repmat(1:4,4,1),1,[]);
  D = sparse(II,JJ,-EV*Q./[A;A;A],4*size(F,1),4*size(V,1));

end

