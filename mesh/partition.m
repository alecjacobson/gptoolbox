function [G,S,D] = partition(W,k,S,WS)
  % PARTITION partition vertices into groups based on each
  % vertex's vector: vertices with similar coordinates (close in 
  % space) will be put in the same group.
  %
  % G = partition(W,k)
  %
  % Inputs:
  %   V  #V by dim coordinate matrix
  %   k  desired number of groups default is dim
  %   S  list of initial seed vertices, guaranteed to be first entries in S
  %   WS  list of initial seed weights, guaranteed to be first entries in S
  % Output:
  %   G  #V list of group indices (1 to k) for each vertex, such that vertex i 
  %     is assigned to group G(i)
  %   S  k  list of seed vertices
  %   D  #V list of squared distances for each vertex to it's corresponding
  %     closest seed
  %
  % Example #1:
  %   % mesh (V,F), weights W
  %   % number of groups
  %   k = 5;
  %   % Compute parititons
  %   [G,S,D] = partition(W,k);
  %   % display mesh colored by partitions
  %   subplot(1,2,1);
  %   trisurf(F,V(:,1),V(:,2),V(:,3),G); 
  %   hold on;
  %   scatter(V(S,1),V(S,2),V(S,3), ...
  %     'o','MarkerFaceColor','y', ...
  %     'MarkerEdgeColor','k',...
  %     'LineWidth',2,'SizeData',200);  
  %   hold off;
  %   view(2);axis equal
  %   title('parition.m');
  %   % Compute parititons using kmeans
  %   [G,Sp,~,D] = kmeans(W,k);
  %   subplot(1,2,2);
  %   % display mesh colored by partitions
  %   trisurf(F,V(:,1),V(:,2),V(:,3),G); 
  %   view(2);axis equal
  %   title('kmeans');
  %
  % Example #2:
  %   % mesh (V,F), weights W
  %   % number of domain vertices
  %   n = size(W,1);
  %   % number of weight functions
  %   m = size(W,2);
  %   % number of groups
  %   k = round(sqrt(n));
  %   % Compute parititons
  %   [G,S,D] = partition(W,k);
  %   % compute averages per group (not the same as seeds)
  %   avgG = ...
  %     sparse(repmat(G,m,1),reshape(repmat(1:m,n,1),n*m,1),W(:))./ ...
  %     repmat(sparse(G,1,1),1,m);
  %   % compute sum of sum of distances per cluster to average
  %   sum(sum((avgG(G,:)-W).^2,2)) 
  %   % use as initial guess for kmeans clustering
  %   [G,S,X] = kmeans(W,k,'Start',W(S,:));sum(X)
  %   % compare to 100 initial guesses
  %   [G,S,X] = kmeans(W,k,'Replicates',100);sum(X)
  %  
  %
  % See also: kmeans
  %



  % number of vertices
  n = size(W,1);
  % number of dimensions 
  m = size(W,2);

  if ~exist('k','var')
    k = m;
  end

  assert(k >= 1);

  % seed vertices
  if ~exist('S','var')
    S = [];
  end

  if ~exist('WS','var')
    WS = W(S,:);
  end

  if isempty(S)
    % "randomly" choose first seed
    % pick a vertex farthest from 0
    [s,I] = sort(sum(W.^2,2),'descend');
    S(1) = I(1);
  end

  % initialize min (squareed) distance to closest seed in seed group
  D = sum((W-repmat(W(S(1),:),n,1)).^2,2);
  % intilialize indices of closest seed
  G = ones(n,1);

  % greedily choose k seeds
  for ii = 2:k
    if ii > numel(S)
      % find maximum
      [s,I] = sort(D,'descend');
      S(ii) = I(1);
      WSii = W(S(ii),:);
    else 
      WSii = WS(ii,:);
    end
    % update min (squared) distance to seed group
    [D,C] = min([D sum((W-repmat(W(S(ii),:),n,1)).^2,2)],[],2);
    % update seed ids for all vertices whose closest seed is the new one
    G(C == 2) = ii;
  end

  [G,centroids,~,D] = kmeans(W,k,'Start',W(S,:));
  % find closest
  [~,S] = min(all_pairs_distances(W,centroids),[],1);
  %warning('not executing kmeans');
  %warning('not executing kmeans');

end
