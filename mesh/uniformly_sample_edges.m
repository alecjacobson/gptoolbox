function S = uniformly_sample_edges(V,E,h)
  % UNIFORMLY_SAMPLE_EDGES Sample edges of (V,E) with spacing close to h.
  %
  % Inputs:
  %   V  #V by dim list of vertex locations
  %   E  #E by 2 list of edges indices into V
  %   h  desired distance between samples
  % Outputs:
  %   S  #S by dim list of sample locations
  %
  S = [];
  for e = 1:size(E,1)
    n = ceil(sqrt(sum((V(E(e,1),:)-V(E(e,2),:)).^2,2))/h)+1;
    s = linspace(0,1,n)';
    S = [S; ...
        bsxfun(@times,1-s,V(E(e,1),:))+bsxfun(@times,  s,V(E(e,2),:))];
  end
  S = unique(S,'rows');
end
