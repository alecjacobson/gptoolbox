function [Imc,A] = apply_matcap(tsh,mc)
% [Imc,A] = apply_matcap(tsh,mc)
%
% Inputs:
%   tsh  handle to a triangle mesh trisurf/patch object
%   mc  mc-height by mc-width by 3 rgb texture image
% Outputs:
%   IO  gcf-height by gcf-width by 3 rgb image of current frame with mesh
%     but face colors are replaced by texture color
%   A  gcf-height by gcf-widsth boolean image of where the mesh (faces) are.
%  
  if ~isfloat(mc)
    mc = im2double(mc);
  end
  face_based = isempty(tsh.VertexNormals);
  if face_based
    N = -tsh.FaceNormals;
  else
    N = -tsh.VertexNormals;
  end
  psh.FaceColor = tsh.FaceColor;
  psh.FaceLighting = tsh.FaceLighting;
  tsh.FaceColor = 'w';
  tsh.FaceLighting = 'none';
  [I,A,FI,B] = shader(tsh);
  [II,~,IV] = find(FI(:));
  V = tsh.Vertices;
  F = tsh.Faces;
  [AZ,EL] = view;
  M = eye(3,4)*viewmtx(AZ,EL)'*eye(4,3);
  %N = per_vertex_normals(V,F)*M;
  N = N*M;
  B1 = B(:,:,1);B2 = B(:,:,2);B3 = B(:,:,3);
  PN = zeros(numel(B1),3);
  if face_based
    PN(II,:) = normalizerow( N(IV,:));
  else
    PN(II,:) = normalizerow( ...
      N(F(IV,1),:).*B1(II) + N(F(IV,2),:).*B2(II) + N(F(IV,3),:).*B3(II));
  end
  PN = reshape(PN,size(B)).*A;
  
  [Xmc,Ymc] = meshgrid(linspace(-1,1,size(mc,2)),linspace(-1,1,size(mc,1)));
  Imc = [];
  Imc(:,:,1) = interp2(Xmc,Ymc,mc(:,:,1),PN(:,:,1),-PN(:,:,2)).*A + I(:,:,1).*~A;
  Imc(:,:,2) = interp2(Xmc,Ymc,mc(:,:,2),PN(:,:,1),-PN(:,:,2)).*A + I(:,:,2).*~A;
  Imc(:,:,3) = interp2(Xmc,Ymc,mc(:,:,3),PN(:,:,1),-PN(:,:,2)).*A + I(:,:,3).*~A;
  tsh.FaceColor = psh.FaceColor;
  tsh.FaceLighting = psh.FaceLighting;
end

