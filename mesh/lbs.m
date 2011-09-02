function [U] = lbs(V,T,W)
  % LBS Compute linear blend skinning deformation of vertices V, using
  % transformations at some control points T, propogated to the mesh using
  % weights W.
  %
  % [U] = lbs(V,T,W)
  % 
  % Inputs:
  %  V  list of vertex positions
  %  T  list of transformations for each controls point, for 2D:
  %    2 by 3 by #controls, for 3D: 3 x 4 by # controls
  %  W  weights, # vertices by # handles matrix of weights
  % Output:
  %  U  list of new vertex positions
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: dualquatlbs
  %

  %% only handle 2D for now
  %assert(size(V,2) == 2);
  dim = size(V,2);

  % number of control points
  m = size(T,3);
  % stack T is one tall 2*m by 3 matrix
  TT = reshape(permute(T,[2,1,3]),[dim+1,dim*m])';
  % Transpose V and append homogeneous w-coordinate
  VV = [V(:,1:dim)'; ones(1,size(V,1))];
  % Multiply each transformation against each mesh vertex position, and reshape
  % so that 3 dim gives deformation according to each control point
  VVV = reshape((TT*VV)',[size(V,1),dim,m]);
  % stack weights for each coordinate for each control
  WW = permute(repmat(W,[1,1,dim]),[1,3,2]);
  % added weighted combinations of deformations according to each control
  U = sum(WW.*VVV,3);

end
