function nim = tangent_to_object_normal_map(V,F,TV,TF,rim)
  % Convert tangent space normal map image to an object space normal map.
  %
  % nim = tangent_to_object_normal_map(V,F,TV,TF,rim);
  %
  % Inputs:
  %   V  #V by 3 list of 3D vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   TV  #TV by 3 list of 2D uv vertex positions
  %   TF  #F by 3 list of triangle indices into TV
  %   rim  #rim by #rim by 3 image of *tangent-space* normals (in range [0,1]),
  %     i.e., a "purple/blue" relative normal map image.
  % Outputs:
  %   nim  #rim by #rim by 3 image of *object-space* normals (in range [0,1])
  %

  [X,Y] = meshgrid(linspace(0,1,size(rim,2)),linspace(1,0,size(rim,1)));
  Q = [X(:) Y(:)];
  I = in_element_aabb(TV(:,1:2),TF,Q);
  % recover barycentric coordinates (assuming tet mesh)
  [BI,~,IV] = find(I);
  W = barycentric_coordinates( ...
    Q(BI,:),TV(TF(IV,1),1:2),TV(TF(IV,2),1:2),TV(TF(IV,3),1:2));
  %TB = sparse(repmat(1:size(W,1),3,1)',TF(IV,:),W,size(W,1),size(TV,1));
  B = sparse(repmat(1:size(W,1),3,1)' , F(IV,:),W,size(W,1),size(V,1));
  R = reshape(rim,[],3).*[2 2 1]-[1 1 0];
  [fT,fN,fB] = per_vertex_frames(V,F,TV,TF);
  NN = normalizerow(sum(cat(3,B*fT,B*fB,B*fN).*permute(R(BI,:),[1 3 2]),3));
  % Use phong normal
  VN = per_vertex_normals(V,F);
  BN = normalizerow(B*VN);
  NN = normalizerow(sum(cat(3,B*fT,B*fB,BN).*permute(R(BI,:),[1 3 2]),3));
  %NN = normalizerow(sum(cat(3,normalizerow(B*fT),normalizerow(B*fB),normalizerow(B*fN)).*permute(R(BI,:),[1 3 2]),3));

  nim = zeros(size(rim,1)*size(rim,2),3);
  nim(BI,:) = NN;
  % Laplacian hole filling / inpainting
  L = fd_laplacian(size(rim(:,:,1)));
  nim = min_quad_with_fixed(-L,[],BI,NN,[],[],struct('force_Aeq_li',true));
  nim = reshape(nim*0.5+0.5,size(rim));
end
