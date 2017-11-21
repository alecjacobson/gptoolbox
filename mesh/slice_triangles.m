function [U,E,J] = slice_triangles(V,F,plane,varargin)
  % SLICE_TRIANGLES 
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangles indices into V
  %   plane  4-long plane equation [nx ny nz -d]
  % Outputs:
  %   U  #U by 3 list of vertex positions
  %   E  #E by 2 list of edge indices into U
  %   J  #E by 1 list of birth triangle indices into F
  % 

  % HACK!
  [U,E,J] = slice_tets(V,F(:,[1 2 3 3]),plane);
  % matlab fails to reutrn one face as row vecto
  if size(E,2) == 1 && size(E,1) == 3
    E = E';
  end
  if isempty(E)
    U = [];
    J = [];
    return
  end
  [U,~,I] = remove_duplicate_vertices(U,eps);
  size(E)
  E(E(:,2) == E(:,1),2) = E(E(:,2) == E(:,1),3);
  E = E(:,1:2);
  N = normals(V,F(J,:));
  N = N-sum(N.*plane(1:3),2).*plane(1:3);
  M = U(E(:,2),:)-U(E(:,1),:);
  Q = [1 plane(1:3)].*[cos(pi/4) sin(pi/4)*[1 1 1]];
  W = quatmultiply(quatmultiply(Q,[zeros(size(M,1),1) M]),Q.*[1 -1 -1 -1]);
  W = W(:,2:4);
  R = sign(sum(W.*N,2))>0;
  E(R,:) = fliplr(E(R,:));
  E = unique(E,'rows');

end
