function [U,G] = unzip_corners(A)
  % UNZIP_CORNERS Given a triangle mesh where corners of each triangle index
  % different matrices of attributes (e.g. read from an OBJ file), unzip the
  % corners into unique efficiently: attributes become properly vertex valued
  % (usually creating greater than #V but less than #F*3 vertices).
  %
  % Inputs:
  %   A  #F by 3 by #A, typically A = cat(3,F,FTC,FN)
  % Outputs:
  %   U  #U by #A list of indices into each attribute for each unique mesh
  %     vertex: U(v,a) is the attribute index of vertex v in attribute a.
  %   G  #F by 3 list of triangle indices into U
  % Example:
  %   [V,F,TC,FTC] = readOBJ('~/Downloads/kiwis/kiwi.obj');
  %   [U,G] = unzip_corners(cat(3,F,FTC));
  %   % display mesh
  %   tsurf(G,V(U(:,1),:));
  %   % display texture coordinates
  %   tsurf(G,TC(U(:,2),:));
  %

  % C  #F*3 by 2 list of vertex indices and uv indices
  C = reshape(A,[],size(A,3));
  [U,~,J] = unique(C,'rows');
  G = J(reshape(1:size(A,1)*3,[],3));
end
