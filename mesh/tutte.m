function U = tutte(V,F)
  % TUTTE Find the longest boundary, fill all other boundaries, and compute a
  % mapping to the plane so that all triangles have positive signed area.
  %
  % U = tutte(V,F) 
  % 
  % Inputs:
  %   V  #V by dim list of input vertex positions
  %   F  #F by 3 list of triangle indices into rows of V
  % Outputs:
  %   U  #V by 2 list of 2D output vertex positions
  %

  % find boundaries
  [O,J] = boundary_faces(F);
  l = edge_lengths(V,O);
  % find longest boundary
  [~,C] = conncomp(facet_adjacency_matrix(O));
  Cl = accumarray(C',l);
  [~,longest] = max(Cl);
  % fill all the other boundaries
  HF = fill_holes(V,F,'Skip',O(find(C==longest,1),1));
  FHF = [F;HF];
  % Only keep longest boundary
  O = O(C==longest,:);
  % vertices in order around boundary
  b = graphtraverse(adjacency_matrix(O),O(1));
  % parametric distance along boundary
  D = cumsum([0;normrow(V(b,:)-V(b([2:end 1]),:))]);
  D = D(1:end-1)/D(end);
  % points on circle to map to
  bc = [cos(D*2*pi) sin(D*2*pi)];
  %% Failed attempt to map to convex hull
  %[~,A] = affine_fit(V(b,:));
  %CV = V(b,:)*A;
  %Cb = convhull(CV);
  %[~,bi,Cbi] = intersect(b,b(Cb));
  %bi = bi(1);
  %Cbi = Cbi(1);
  %b = b([bi:end 1:bi-1]);
  %D = D([bi:end 1:bi-1]);
  %Cb = Cb([Cbi:end-1 1:Cbi-1 Cbi]);
  %CD = cumsum([0;normrow(CV(Cb,:)-CV(Cb([2:end 1]),:))]);
  %CD = CD(1:end-1)/CD(end);
  %CV = CV(Cb,:);
  %D = max(D)-D;
  %[~,K] = histc(D,CD);
  %T = (D-CD(K))./(CD(K+1)-CD(K));
  %bc = CV(K,:) + T.*(CV(K+1,:)-CV(K,:));
  
  % scale so that area matches
  total_area = sum(doublearea(V,F))*0.5;
  bc = bc*sqrt(total_area/pi);
  % Cascading attempts to build mapping to disk
  % cotangent Laplacian
  U = kharmonic(V,FHF,b,bc,1);
  if min(doublearea(U,F)) < 0
    % Intrinsic Laplacian. This "should" work, but due to numerics, might not.
    U = kharmonic(V,FHF,b,bc,1,'IntrinsicDelaunay',true);
  end
  if min(doublearea(U,F)) < 0
    U = kharmonic([],FHF,b,bc,1);
  end
  assert(min(doublearea(U,F))>0);
end
