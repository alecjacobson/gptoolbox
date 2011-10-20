function [near,in,on,B,L] = near_mesh(V,F,Q,epsilon)
  % NEAR_MESH test whether a list of points are in or near a given mesh
  %  
  % [near] = near_mesh(V,F,Q)
  % [near,in,on,B,L] = near_mesh(V,F,Q)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of face indices
  %   Q  #Q by dim list of query points
  %   epsilon  minium distance allowed after collapses are complete, default is to
  %     use fraction of maximum edge length
  % Outputs:
  %   near #Q list of flags revealing whether queries are near (V,F)
  %   in #Q list of flags revealing whether queries are in (V,F)
  %   on #Q list of flags revealing whether queries are on boundary of (V,F)
  %   B  #B by 1 list of mesh outline edges 
  %   L  #loops+1 by 1 list of boundary loop start indices into B, the last
  %     entries is (by tradition) always the numel of B + 1
  %
  % See in_mesh, inpolygon
  %

  if ~exist('epsilon','var') || isempty(epsilon)
    EE = edges(F);
    % maximum edge length
    %maxD = max(sqrt(sum((V(EE(:,1),:) - V(EE(:,2),:)).^2,2)));
    minD = min(sqrt(sum((V(EE(:,1),:) - V(EE(:,2),:)).^2,2)));
    epsilon = minD/2;
  end

  dim = size(V,2);
  % only works in 2D
  assert(dim == 2);

  % first determine points strictly in or on mesh
  [in,on,B,L] = in_mesh(V,F,Q);

  % avoid sqrts
  sqr_eps = epsilon.^2;

  % boundary edges
  BE = [B; B(2:end) B(1)]';
  % compute projection of each point to each boundary line segment
  [T,sqrD] = project_to_lines(Q,V(BE(:,1),:),V(BE(:,2),:));
  % each vertex seen by each edge
  QBE = repmat(Q,[1 1 size(BE,1)]);
  % edge start positions
  S = V(BE(:,1),:);
  % edge destination positions
  D = V(BE(:,2),:);
  % distance of each point to each edge start
  sqrDS = ...
    squeeze(sum((QBE - permute(repmat(S,[1 1 size(Q,1)]),[3 2 1])).^2,2));
  % distance of each point to each edge dest
  sqrDD = ...
    squeeze(sum((QBE - permute(repmat(D,[1 1 size(Q,1)]),[3 2 1])).^2,2));
  % replace distances to edges when point is closest to start or dest endpoints
  % respectively
  sqrD(T<0) = sqrDS(T<0);
  sqrD(T>1) = sqrDD(T>1);
  % compute minimum distance to boundary
  [minD] = min(sqrD,[],2);
  % mask telling whether closest point for each edge is close enough
  near = minD<sqr_eps;
  % all strictly in points are also close
  near = in | near;
end
