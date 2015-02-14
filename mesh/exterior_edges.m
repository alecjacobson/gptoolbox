function [E,C,aC] = exterior_edges(F)
  % EXTERIOR_EDGES List of exterior (boundary of manifold patch) edges. Note:
  % this is---in general---not the same as the _open boundary_ of the mesh.
  %
  % Inputs:
  %   F  #F by dim=3 list of facet indices
  % Outputs:
  %   E  #E by 2 list of exterior edges
  %   C  #E by 1 list of signed Counts
  %   aC  #E by 1 list of unsigned Counts
  %
  % See also: outline
  %

  switch size(F,2)
  case 2
    E = F;
    % TODO: This is not robust to non-manifold issues
    % edge going out
    out =  sparse(E(:,1),1,1,max(E(:)),1);
    % edge coming in
    in  = sparse(E(:,2),1,1,max(E(:)),1);
    % endpoints
    ep = (in + out) == 1;
    epi = find(ep);
    E = [repmat(epi(1),numel(epi)-1,1) epi(2:end)];
    E = ...
      bsxfun(@times,out(E(:,2))   ,E) +  ...
      bsxfun(@times,(~out(E(:,2))),fliplr(E));
  case 3
    allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
    sortallE = sort(allE,2);
    % -1 in this direction +1 in that direction --> stored in this direction
    C = sparse(sortallE(:,1),sortallE(:,2),(allE(:,1)<allE(:,2))*2-1);
    % +1 in any direction --> stored in this direction
    aC = sparse(sortallE(:,1),sortallE(:,2),1);
    assert(nnz(tril(C)) == 0);
    [~,~,aC] = find((C~=0).*aC);
    [EI,EJ,C] = find(C);
    E = [EI EJ];
  end

end
