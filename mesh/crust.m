function CE = crust(P)
  % CRUST Given a set of points P determine the set of facets CE that form the
  % "crust" (assuming points are a sampling of some underlying curve/surface).
  % Implements "A New Voronoi-Based Surface Reconstruction Algorithm" [Amenta
  % et al. 1998].
  %
  % Inputs:
  %   P  #P by dim list of input points
  % Outputs:
  %   CE  #CE by dim list of facets indexing P
  %
  dim = size(P,2);
  switch dim
  case 2
    % As of at least MATLAB 2017a, delaunayTriangulation is actually faster
    % than calling triangle
    %
    %[~,F] = triangle(P,[],[],'Flags','-c');
    Tri = delaunayTriangulation(P);
    F = Tri.ConnectivityList;
    [~,C] = circumradius(P,F);
    C = C(~any(isnan(C),2),:);
    C = unique(C,'rows');
    PC = [P;C];
    %[~,FF] = triangle(PC,[],[],'Flags','-c');
    Tri = delaunayTriangulation(PC);
    FF = Tri.ConnectivityList;
    E = edges(FF);
    CE = E(all(E<=size(P,1),2),:);
  end
end
