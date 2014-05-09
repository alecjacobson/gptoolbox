function D = all_pairs_distances(V,U)
  % ALL_PAIRS_DISTANCES compute distances between each point i in V and point j
  % in U
  % 
  % This is obsolete: use pdist2 instead
  % 
  % D = all_pairs_distances(V,U)
  % 
  % Inputs:
  %   V  #V by dim list of points
  %   U  #U by dim list of points
  % Outputs:
  %   D  #V by #U matrix of distances, where D(i,j) gives the distance between
  %     V(i,:) and U(j,:)
  % 

  warning('obsolete. Call `pdist2` directly instead');

  assert(size(V,2) == size(U,2));

  %% Super-slow for loops method
  %% Elapsed time is 128.969203 seconds.
  %D = zeros(size(V,1),size(U,1));
  %for i = 1:size(V,1)
  %  for j = 1:size(U,1)
  %    D(i,j) = sqrt(sum((V(i,:)-U(j,:)).^2,2));
  %  end
  %end

  %% Slow, especially if hitting memory paging
  %% Elapsed time is 2.519262 seconds.
  %D = ...
  %  squeeze(sqrt(sum((repmat(V,[1 1 size(U,1)])- ...
  %  permute(repmat(U,[1 1 size(V,1)]),[3 2 1])).^2,2)));

  %% Faster single for loop, single repmat method, handles memory paging better
  %% Elapsed time is 1.160419 seconds.
  %D = zeros(size(V,1),size(U,1));
  %for i = 1:size(V,1)
  %  D(i,:) = ...
  %    sqrt(sum((repmat(V(i,:),[size(U,1) 1]) - U).^2,2))';
  %end



  if size(V,1) > 1000 || size(U,1) > 1000
    % Fast single for loop, single bsxfun method, but matlab can handle both for
    % loops in bsxfun. This seems to be a nice balance between hitting a memory
    % wall for big input and still being fast for medium and slow input
    % Elapsed time is 0.863138 seconds.
    D = zeros(size(V,1),size(U,1));
    for i = 1:size(V,1)
      D(i,:) = ...
        sqrt(sum(bsxfun(@minus,V(i,:),U).^2,2));
    end
  else
    % Fastest
    % Elapsed time is 0.653082 seconds.
    D = permute(sqrt(sum(bsxfun(@minus,V,permute(U,[3 2 1])).^2,2)),[1 3 2]);
  end
end
