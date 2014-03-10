function [B,L] = ordered_outline(F)
  % ORDERED_OUTLINE Find outline (boundary) edges of mesh
  %
  % [B,L] = outline(F)
  %
  % Input:
  %   F  #F by 3 face list of indices
  % Outputs:
  %   B  #B by 1 list of outline edges 
  %   L  #loops+1 by 1 list of boundary loop start indices into B, the last
  %     entries is (by tradition) always the numel of B + 1
  %
  % Example:
  %   V = [0 0; 1 0; 1 1 ; 0 1; 4 0; 4 4; 0 4];
  %   F = [1 2 3; 1 2 4; 6 5 7];
  %   [B,L] = ordered_outline(F);
  %   hold on
  %   for l = 1:(numel(L)-1)
  %     plot(V([B(L(l):(L(l+1)-1)) B(L(l))],1),V([B(L(l):(L(l+1)-1)) B(L(l))],2));
  %   end
  %   hold off
  %
  % See also: outline
  %

  O = outline(F);
  % determine uniqueness of indices 
  [u,m,n] = unique(O(:),'rows');
  % determine counts for each unique edge
  counts = accumarray(n(:), 1);
  % all counts should be 2
  assert(all(counts == 2));

  % number of outline edges
  no = size(O,1);

  % list of boundary 
  B = [];
  L = [];
  unseen = true(no,1);
  while any(unseen)
    L = [L numel(B)+1];
    % find first unseen index
    f = find(unseen,1);
    % append start vertex to boundary
    next = O(f,1);
    % keep track of start
    start = next;
    unseen(f) = false;
    B = [B next];
    next = O(f,2);
    % loop until back at start
    while next ~= start
      B = [B next];
      % try to find next as source
      f = find(unseen & O(:,1) == next);
      if isempty(f)
        % try to find next as dest
        f = find(unseen & O(:,2) == next);
        assert(~isempty(f));
        next = O(f,1);
      else
        next = O(f,2);
      end
      unseen(f) = false;
    end
  end
  L = [L numel(B)+1];

end

