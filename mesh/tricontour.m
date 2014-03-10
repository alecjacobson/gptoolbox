function tricontour(V,F,w,speedup)
  % TRICONTOUR Plot a filled contour defined on top of a triangle mesh returns
  % function handle
  %
  % tricontour(V,F,w,speedup)
  %
  %
  % Inputs:
  %  V  list of vertex positions
  %  F  list of face indices
  %  w  scalar values defined on V
  %  speedup optional speed up by not eliminating vertices outside shape {false}

  if(~exist('speedup','var'))
    speedup = false;
  end
  [Xr,Yr,Wr] = triinterp(V,F,w,speedup);
  contourf(Xr,Yr,Wr)

  % Find all edges in mesh, note internal edges are repeated
  E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
  % determine uniqueness of edges
  [u,m,n] = unique(E,'rows');
  % determine counts for each unique edge
  counts = accumarray(n(:), 1);
  % extract edges that only occurred once
  O = u(counts==1);

  line(V([O(:);O(1)],1),V([O(:);O(1)],2),'Color','k','LineWidth',6)
end
