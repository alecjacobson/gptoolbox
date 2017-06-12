function [CE,PC] = crust(P)
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
  case 3
    % As of MATLAB 2017a, much faster than calling tetgen
    Tri = delaunayTriangulation(P);
    [VV,VR] = Tri.voronoiDiagram;
    n = size(P,1);
    I = zeros(n,2);
    %% This loop is the bottleneck ~50% but actually scales reasonable well with
    %% input size (as expected since it is O(n) and the delaunay is O(n log n) )
    %for vi = 1:n
    %  VRi = VR{vi};
    %  VVi = VV(VRi,:);
    %  D = pdist2(VVi,P(vi,:));
    %  Ii = [0 0];
    %  [~,Ii(1)] = max(D);
    %  D = D.*sign(sum((VVi(Ii(1),:)-P(vi,:)).*(VVi-P(vi,:)),2));
    %  [~,Ii(2)] = min(D);
    %  I(vi,:) = VRi(Ii);
    %end

    % Vectorized version of for loop
    counts = cellfun(@(r) numel(r),VR);
    % Creepy way of getting indices repeated variable number of times
    J = cumsum([1;sparse(cumsum(counts),1,1)])';
    J = J(1:end-1);
    VRA = [VR{:}];
    U = VV(VRA,:) - P(J,:);
    D = normrow(U);
    S = sparse(VRA,J,D,size(VV,1),size(P,1));
    [~,I(:,1)] = max(S);
    Y = VV(I(:,1),:)-P;
    D = sign(sum(U.*Y(J,:),2)).*D;
    S = sparse(VRA,J,D,size(VV,1),size(P,1));
    [~,I(:,2)] = minnz(S);

    I = unique(I(:));
    C = VV(I,:);
    C = C(~any(isnan(C)|isinf(C),2),:);
    PC = [P;C];
    Tri = delaunayTriangulation(PC);
    FF = Tri.ConnectivityList;
    E = [ ...
      FF(:,[4 2 3]); ...
      FF(:,[3 1 4]); ...
      FF(:,[2 4 1]); ...
      FF(:,[1 3 2])];
    E = unique(sort(E,2),'rows');
    CE = E(all(E<=size(P,1),2),:);
  case 2
    % As of at least MATLAB 2017a, delaunayTriangulation is actually faster
    % than calling triangle
    %
    %[~,F] = triangle(P,[],[],'Flags','-c');
    %F = Tri.ConnectivityList;
    %[~,C] = circumradius(P,F);
    %C = C(~any(isnan(C),2),:);
    %C = unique(C,'rows');
    Tri = delaunayTriangulation(P);
    C = Tri.circumcenter;
    PC = [P;C];
    %[~,FF] = triangle(PC,[],[],'Flags','-c');
    Tri = delaunayTriangulation(PC);
    FF = Tri.ConnectivityList;
    E = edges(FF);
    CE = E(all(E<=size(P,1),2),:);
  end
end
