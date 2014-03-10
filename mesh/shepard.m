function [W,dist_H_V] = shepard(V,C,P,E,CE,p)
  % SHEPARD Compute shepard weights for a list of vertices, given a list of
  % samples and optionally a denominator power value.
  %
  % W = shepard(V,C)
  %
  % Inputs:
  %  V  list of vertex positions
  %  C  list of control vertices
  %  P  list of indices into C for point controls, { 1:size(C,1) }
  %  E  list of bones, pairs of indices into C, connecting control vertices, 
  %    { [] }
  %  CE  list of "cage edges", pairs of indices into ***P***, connecting
  %    control ***points***. A "cage edge" just tells point boundary conditions 
  %    to vary linearly along straight lines between the end points, and to be
  %    zero for all other handles. { [] }
  %  p  (optional) power for denominator, scalar or list same size as C {2}
  %
  % Outputs:
  %  W  weights, # vertices by # handles matrix of weights
  %
  % Example:
  %  % "Shape-aware" shepard weights
  %  D = permute(sum((repmat(V,[1,1,size(C,1)]) - ...
  %    permute(repmat(C,[1,1,n]),[3,2,1])).^2,2),[1,3,2]);
  %  [minD,b] = min(D);
  %  snap_C = V(b,:);
  %  [BV,ev,ed] = biharmonic_embedding(V,F,6);
  %  [W] = shepard(BV,BV(b,:));
  %
  %

  % check if either are 3D but really all z's are 0
  V_flat = size(V,2) == 3 && (sqrt(sum(V(:,3).^2)) < 1e-10);
  C_flat = size(C,2) == 3 && (sqrt(sum(C(:,3).^2)) < 1e-10);
  % is both are essentially 2D then ignore z-coords
  if((size(C,2) == 2 || C_flat) && (size(V,2) == 2 || V_flat))
    % ignore z coordinate
    V = V(:,1:2);
    C = C(:,1:2);
  end

  assert(size(C,2) == size(V,2));
  dim = size(C,2);

  % default p value
  if(~exist('p','var'))
    p = 2;
  end

  if(~exist('P','var'))
    P = 1:size(C,1);
  end

  if(~exist('E','var'))
    E = [];
  end

  if(~exist('CE','var'))
    CE = [];
  end

  assert(prod(size(P)) == max(size(P)));
  assert(isempty(E) || size(E,2) == 2);
  assert(isempty(CE) || size(CE,2) == 2);
  assert(isempty(E) || max(E(:))<=size(C,1));
  assert(isempty(CE) || max(CE(:))<=numel(P));
  assert(isempty(P) || max(P(:))<=size(C,1));
  assert(isempty(E) || min(E(:)) >= 1);
  assert(isempty(CE) || min(CE(:)) >= 1);
  assert(isempty(P) || min(P(:)) >= 1);


  % number of point controls 
  np = numel(P);
  nb = size(E,1);
  % number of domain vertices
  n = size(V,1);

  % if p is a scalar convert it now to a list the same size as C
  if(prod(size(p)) == 1)
    p = repmat(p,np+nb,1);
  elseif ~isempty(CE)
    if 0 ~= std(p(CE(:)))
      error('Variable powers not supported with cage edges');
    end
  end

  dist_P_V = zeros(np,n);
  if(np > 0)
    % vectors from V to every P, where PmV(i,j,:) is the vector from domain
    % vertex j to handle i
    PmV = ...
      permute( ...
        permute(repmat(C(P,:),[1,1,n]),[3,2,1]) - ...
        repmat(V,[1,1,np]),[3,1,2]);
    % distance from V to every P, where dist_P_V(i,j) is the distance from domain
    % vertex j to point handle i
    dist_P_V = sqrt(sum(PmV.^2,3));
    % distance from each corner in P to the next corner so that edge_length(i) 
  end

  dist_B_V = zeros(nb,n);
  if(nb > 0)
    % loop over bones
    for( ii = 1:nb )
      % project to line
      [t,sqr_d] = project_to_lines(V,C(E(ii,1),:),C(E(ii,2),:));
      d = sqrt(sqr_d);
      % if projected point falls outside of line segment take distance to
      % closest end point
      d(t<0) = sqrt(sum((V(t<0,:)-repmat(C(E(ii,1),:),size(V(t<0,:),1),1)).^2,2));
      d(t>1) = sqrt(sum((V(t>1,:)-repmat(C(E(ii,2),:),size(V(t>1,:),1),1)).^2,2));
      dist_B_V(ii,:) = d';
    end
  end

  dist_H_V = [dist_P_V ; dist_B_V];
  % power of each control point seen by each vertex in domain
  pp = repmat(p,1,n);
  W = 1.0./((dist_H_V).^pp);

  if(size(CE,1) > 0)
    % zero out previous entries in points with incident cage edges
    W(CE(:),:) = 0;
    % compute projection of each point to each line segment
    [T,sqrD] = project_to_lines(V,C(P(CE(:,1)),:),C(P(CE(:,2)),:));
    % compute weights for each cage edge as if it were a bone
    pce = (p(CE(:,1))+p(CE(:,2)))/2;
    [~,dist_CE_V] = shepard(V,C,[],P(CE),[],pce);
    % compute "bone weights" for each cage edge as they were before
    % normalization
    ppce = repmat(pce,1,n);
    WCE = 1.0./((dist_CE_V).^ppce)';
    % clamp distances
    T(T<0) = 0;
    T(T>1) = 1;
    % multiply "bone weight" by projection parameter and sum up
    II = repmat(1:n,1,2*size(CE,1))';
    JJ = reshape(repmat(CE(:)',n,1),[],1);
    VV = [WCE(:).*(1-T(:)); WCE(:).*T(:)];
    W = W+sparse(II,JJ,VV,n,np+nb)';
    % cage edges on-samples
    on_sample = dist_CE_V < eps;
    [I,J] = find(on_sample);
    W(:,any(on_sample)) = 0;
    % grab T transpose
    TT = T';
    W(sub2ind(size(W),CE(I,1),J)) = 1-TT(on_sample);
    W(sub2ind(size(W),CE(I,2),J)) = TT(on_sample);
  end

  % Handle degenrate case that a control point is on a mesh vertex
  % snap vertices close to corners
  on_sample = dist_H_V < eps;
  W(:,any(on_sample,1)) = 0;
  W(on_sample) = 1;

  % normalize W
  W = W./repmat(sum(W,1),np+nb,1);

  % we've made W transpose
  W = W';

end
