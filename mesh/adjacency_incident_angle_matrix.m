function [A,AE] = adjacency_incident_angle_matrix(V,E)
  % ADJACENCY_INCIDENT_ANGLE_MATRIX  Compute an adjacency matrix between *edges*
  % where non-zeros are the angles between adjacent edges (edges that share a
  % vertex).
  %
  % [A] = adjacency_incident_angle_matrix(V,E)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   E  #E by 2 list of edge indices into rows of V
  % Outputs:
  %   A  #E by #E adjacency matrix so that A(i,j) = h ≠ 0 implies that edges
  %     E(i,:) and E(j,:) meet at a coincident vertex with angle h.
  %   AE  #E by #E adjacency_matrix so that A(i,j) ≠ 0 implies that edges 
  %     E(i,:) and E(j,:) meet at a coincident vertex (just in case h==0)
  %

  assert(size(E,1)==size(unique(sort(E,2),'rows'),1));

  V2E = sparse(E,repmat(1:size(E,1),2,1)',1,size(V,1),size(E,1));
  AE = V2E'*V2E;
  [I,J] = find(triu(AE,1));
  C = sparse(repmat(1:numel(I),4,1)',[E(I,:) E(J,:)],1,numel(I),size(V,1));
  MI = (1:numel(I))';
  % get center vertex
  [~,M2] = max(C==2,[],2);
  C(sub2ind(size(C),MI,M2)) = 0;
  % get other vertex (who cares which order)
  [~,M1] = max(C,[],2);
  C(sub2ind(size(C),MI,M1)) = 0;
  % get other vertex (who cares which order)
  [~,M3] = max(C,[],2);
  V12 = normalizerow(V(M1,:)-V(M2,:));
  V23 = normalizerow(V(M3,:)-V(M2,:));
  H = acos(sum(V12.*V23,2));
  % rebuild A
  A = sparse([I J],[J I],[H H],size(E,1),size(E,1));
end
