function [SV,SF,VV,FF] = linear_sweep(V,F,sweep)
  % LINEAR_SWEEP Compute the surface of the solid sweep of a surface mesh (V,F)
  % along a vector (sweep): i.e. the Minkowski sum of (V,F) and the line
  % segment (0,0)-->(sweep) 
  %
  % [SV,SF] = linear_sweep(V,F,sweep)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by dim list of facet indices into V
  %   sweep  1 by dim vector defining sweep direction and distance
  % Outputs:
  %   SV  #SV by dim list of vertex positions
  %   SF  #F by dim list of facet indices into SV
  %

  dim = size(V,2);
  assert(dim == size(F,2),'Facet degree should equal dim');

  % vertices at end of sweep
  U = bsxfun(@plus,V,sweep);

  % number of vertices
  n = size(V,1);

  switch dim
  case 2
    % Facets are edges
    E = F;
    % normals
    N = normalizerow(V(E(:,2),:)-V(E(:,1),:))*[0 -1;1 0];
    % dot product with sweep direction
    D = sum(bsxfun(@times,N,sweep),2);
    % edges with positive and negative dot product with sweep
    Ep = E(D>0,:);
    Em = E(D<=0,:);
    % sign of incidence for all vertices
    S = sparse(Ep,1,repmat([1 -1],size(Ep,1),1),n,1);
    % find any incidence
    I = find(S);
    % "walls"
    Epw = [I n+I];
    Epw(S(I)<0,:) = fliplr(Epw(S(I)<0,:));
    EE = [E;n+E;Epw];
    VV = [V;U];
    % Tesselated interior
    [TV,TF] = triangle(VV,EE,[],'Quiet');
    % (self-intersecting) boundary Sweep
    WE = [Em;n+Ep;Epw];
    WV = [V;U];
    % Classify interior elements with winding number
    w = winding_number(WV,WE,barycenter(TV,TF));
    % Extract boundary of w~=0 part
    SE = outline(TF(round(w)~=0,:));
    [SV,IM] = remove_unreferenced(TV,SE);
    SE = IM(SE);
    % output names
    SF =SE;
  case 3
    % normals
    N = normalizerow(normals(V,F));
    % dot product with sweep direction
    D = sum(bsxfun(@times,N,sweep),2);
    % faces with positive and negative dot product with sweep
    Fp = F(D>0,:);
    Fm = F(D<=0,:);
    O = outline(Fp);
    Fpw = [O(:,[1 2]) n+O(:,2);n+O(:,[2 1]) O(:,1)];

    VV = [V;U];
    FFp = [fliplr(Fp);Fpw;n+Fp];
    FFm = [Fm;Fpw;n+fliplr(Fm)];

    [VVp,IM] = remove_unreferenced(VV,FFp);
    FFp = IM(FFp);
    [VVm,IM] = remove_unreferenced(VV,FFm);
    FFm = IM(FFm);
    writePLY('front.ply',VVp,FFp);
    writePLY('back.ply',VVm,FFm);
    statistics(VVp,FFp)
    statistics(VVm,FFm)
    [SV,SF] = mesh_boolean(VVp,FFp,VVm,FFm,'union');

    %% (self-intersecting) boundary Sweep
    %FF = [Fm;n+Fp;Fpw];
    %VV = [V;U];

    %%% Self-union: this doesn't work... Non-manifold edges in Fpw? Yes, this is
    %%% possible (folded Y case), but why is it a problem?
    %%warning('experimental. not working.');
    %%[SV,SF] = mesh_boolean(VV,FF,[],[],'union');
    %%return;

    %% Tesselated interior
    %[SVV,SFF,~,~,IM] = selfintersect(VV,FF);
    %SFF = IM(SFF);
    %[SVV,IM] = remove_unreferenced(SVV,SFF);
    %SFF = IM(SFF);

    %% Can do this all using "peeling" like booleans. For now just do single
    %% outer-layer in easy case.

    %% check if bounding boxes overlap
    %[BV,BF] = bounding_box(V);
    %BU = bsxfun(@plus,BV,sweep);
    %if false && all(round(winding_number(BV,BF,BU))==0)
    %  SF = outer_hull(SVV,SFF);
    %  [SV,IM] = remove_unreferenced(SVV,SF);
    %  SF = IM(SF);
    %else
    %  %statistics(SVV,SFF)
    %  [TV,TT,TF] = tetgen(SVV,SFF,'Flags','-YT1e-16');
    %  BC = barycenter(TV,TT);
    %  % Classify interior elements with winding number
    %  w = winding_number(VV,FF,BC);
    %  numel(w(abs(w-round(w))>1e-7))
    %  % Extract boundary of w~=0 part
    %  SF = boundary_faces(TT(round(w)~=0,:));
    %  [SV,IM] = remove_unreferenced(TV,SF);
    %  SF = IM(SF);
    %end
  end

end
