function [I,J,K] = edges_to_path(E)
  % EDGES_TO_PATH Given a set of undirected, unique edges such that all form a
  % single connected compoent with exactly 0 or 2 nodes with valence =1,
  % determine the/a path visiting all nodes.
  %
  % [I,J,K] = edges_to_path(E)
  %
  % Inputs:
  %   E  #E by 2 list of undirected edges
  % Outputs:
  %   I  #E+1 list of nodes in order tracing the chain (loop), if the output
  %     is a loop then I(1) == I(end)
  %   J  #I-1 list of indices into E of edges tracing I
  %   K  #I-1 list of indices into columns of E {1,2} so that K(i) means that
  %     E(i,K(i)) comes before the other (i.e., E(i,3-K(i)) ). This means that 
  %     I(i) == E(J(i),K(i)) for i<#I, or
  %     I == E(sub2ind(size(E),J([1:end end]),[K;3-K(end)]))))
  % 
  
  if size(E,1) == 1
    I = E(:);
    J = 1;
    K = 1;
    return;
  end

  % Compute on reduced graph
  E_orig = E;
  [U,~,E] = unique(E_orig(:));
  E = reshape(E,size(E_orig));

  % Valences
  V = sparse(E(:),1,1);
  assert(all(V<=2));
  % Try to find a vertex with valence = 1
  [c,s] = minnz(V);
  %I = graphtraverse(adjacency_matrix(E),s)';
  I = reshape(dfsearch(graph(adjacency_matrix(E)),s),[],1);
  if c == 2
    I = I([1:end 1]);
  end

  if nargout>1
    % This is slow. On way to fix this is to store the index into E of the edge
    % connecting node i to node j in the "adjacency matrix" above. Then can
    % just read off values (might be faster)
    [F,J] = ismember(sort([I(1:end-1) I(2:end)],2),sort(E,2),'rows');
    assert(all(F));
  end
  if nargout>2
    K = (E(J,1) ~= I(1:end-1))+1;
    assert(all(I == E(sub2ind(size(E),J([1:end end]),[K;3-K(end)]))));
  end

  % Map vertex indices onto original graph
  I = U(I);
end
