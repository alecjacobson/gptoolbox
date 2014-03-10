function [Xr,Yr,Wr] = triinterp(V,F,w,speedup)
  % TRIINTERP Given a scalar field S defined over a triangle mesh with vertices
  % V and faces F, linearly interpolate the data for new positions U. NaNs are
  % returned if position in U is not in or on a mesh triangle
  %
  % [Xr,Yr,Wr] = triinterp(V,F,w,speedup)
  %
  % Inputs:
  %   V vertex position #V x 3 or #V x 2
  %   F list of faces #F x 3
  %   w scalar field defined over V, #V x 1
  %  speedup optional speed up by not eliminating vertices outside shape {false}
  %
  % Outputs:
  %   Xr  x values for griddata
  %   Yr  y values for griddata
  %   Wr  scalar values for griddata
  %

  [Xr,Yr,Wr] = griddata(V(:,1),V(:,2),w,unique(V(:,1)),unique(V(:,2))');
  if(~exist('speedup','var') || ~speedup)
    % Find all edges in mesh, note internal edges are repeated
    E = sort([F(:,1) F(:,2); F(:,2) F(:,3); F(:,3) F(:,1)]')';
    % determine uniqueness of edges
    [u,m,n] = unique(E,'rows');
    % determine counts for each unique edge
    counts = accumarray(n(:), 1);
    % extract edges that only occurred once
    O = u(counts==1);
    % this is very slow
    IN = inpolygon(Xr,Yr,V(O,1),V(O,2));
    % don't draw points outside the mesh
    Wr(~IN) = NaN;
  end
end
