% IN_ELEMENT_AABB Use an axis-aligned bounding box tree to determine if a set
% of points appears inside the elements of a mesh.
%
% I = in_element_aabb(V,Ele,Q);
% [I,bb_mins,bb_maxs,elements] = ...
%   in_element_aabb(V,Ele,Q,bb_mins,bb_maxs,elements);
%
% Inputs:
%   V  #V by dim list of mesh vertex positions. 
%   Ele  #Ele by dim+1 list of mesh indices into #V. 
%   Q  #Q by dim list of query point positions
%   Optional:
%     bb_mins  max_tree by dim list of bounding box min corner positions
%     bb_maxs  max_tree by dim list of bounding box max corner positions
%     elements  max_tree list of element or (not leaf id) indices into Ele
% Outputs:
%   I  #Q list of indices into Ele of first containing element (0 means no
%     containing element)
%   bb_mins (see optional input)
%   bb_maxs (see optional input)
%   elements (see optional input)
%
% Example:
%   % Dummy data
%   aabbn = [];
%   aabbx = [];
%   aabbe = [];
%   % Call once to build AABB and query
%   [I,aabbn,aabbx,aabbe] = in_element_aabb(V,Ele,Q,aabbn,aabbx,aabbe);
%   % Call again with same output to build AABB from serialized output
%   [I,aabbn,aabbx,aabbe] = in_element_aabb(V,Ele,Q,aabbn,aabbx,aabbe);
%   % recover barycentric coordinates (assuming tet mesh)
%   [II,~,IV] = find(I);
%   B = barycentric_coordinates( ...
%     Q(II,:),V(Ele(IV,1),:),V(Ele(IV,2),:),V(Ele(IV,3),:),V(Ele(IV,4),:));
%   % reproduce positions of found points using barycentric coordinates
%   BQ = sum(bsxfun(@times, permute(B,[1 3 2]), ...
%     cat(3,V(Ele(IV,1),:),V(Ele(IV,2),:),V(Ele(IV,3),:),V(Ele(IV,4),:))),3);
%   % Map back so that size(B,1) == size(Q,1)
%   B = sparse(repmat(II,1,size(B,2)),repmat(1:size(B,2),size(B,1),1),B, ...
%     numel(I),size(B,2));
%   scatter3(Q(:,1),Q(:,2),Q(:,3));
%   hold on;
%   scatter3(BQ(:,1),BQ(:,2),BQ(:,3),'r','SizeData',10);
%   hold off;
%
