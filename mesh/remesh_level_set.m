function [U,G,I,S,LE,H] = remesh_level_set(V,F,D)
  % [U,G,I,S] = remesh_level_set(V,F,D)
  % 
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   D  #V by 1 list of level set values
  % Outputs:
  %   U  #U by dim list of vertex positions
  %   G  #G by 3 list of triangle indices into U
  %   I  #G by 1 list of indices into F
  %   S  #U by #V sparse matrix so that U = S*V
  %   LE  #LE by 2 list of edge indices into rows of U tracing the level set,
  %   oriented CCW around negative regions.
  %  
  % Example:
  %   C = graph_coloring(facet_adjacency_matrix(TF),9);
  %   [U,G,I,S] = remesh_level_set(TV,TF,D);
  %   colormap(cbrewer('Set1',max(C)))
  %   tsurf(G,S*TV,'CData',C(I));

  function local_F = unique_faces(local_T)
    local_allF = [ ...
      local_T(:,[4 2 3]); ...
      local_T(:,[3 1 4]); ...
      local_T(:,[2 4 1]); ...
      local_T(:,[1 3 2])];
    local_sortF = sort(local_allF,2);
    [local_F,~,local_FMAP] = unique(local_sortF,'rows');
  end

  switch size(F,2)
  case 4
    T = F;
    allF = [ ...
      T(:,[4 2 3]); ...
      T(:,[3 1 4]); ...
      T(:,[2 4 1]); ...
      T(:,[1 3 2])];
    sortF = sort(allF,2);
    MF = reshape( ...
      all(allF(:,[1 2 3])==sortF,2) | ...
      all(allF(:,[2 3 1])==sortF,2) | ...
      all(allF(:,[3 1 2])==sortF,2), size(T));
    [F,~,FMAP] = unique(sortF,'rows');
    FMAP = reshape(FMAP,size(T));
    [U,G,I,S,LE,H] = remesh_level_set(V,F,D);


    % J(f,:) = [i j k] means that the f-th face of F was split into triangle i
    % and quad (j,k) further split into triangles j and k in G.
    J = full_sparse(I(H>0),H(H>0),find(H>0),size(F,1),3);

    % go through the motions of trusting recursive output to tel us which tets
    % were not touched at all (rather than just looking at D)
    F_was_split = false(size(F,1),1);
    F_was_split(I) = H>0;
    %T_was_split = false(size(T,1),1);
    TF_was_split = reshape(F_was_split(FMAP),size(T));
    assert(all(sum(TF_was_split,2) == 0 | sum(TF_was_split,2) == 3 | sum(TF_was_split,2) == 4));
    T0 = sum(TF_was_split,2) == 0;
    T3 = sum(TF_was_split,2) == 3;
    T4 = sum(TF_was_split,2) == 4;

    TU = U;
    TG0 = T(T0,:);

    %%B = [];
    %TG3 = [];
    %for i = 1:size(T,1)
    %  if ~T3(i)
    %    continue;
    %  end
    %  % index of vertex on the "cap"
    %  [~,j] = min(TF_was_split(i,:),[],2);
    %  T(i,j)

    %  %clf;
    %  %%tsurf(F,V,'CData',D,falpha(0.5,1),fphong);
    %  %hold on;
    %  %tsurf(edges(G),U,'VertexIndices',1,'CData',S*D,falpha(0.1,1),fphong);
    %  %tsurf(unique_faces(T(i,:)),U,'FaceColor',blue,falpha(0.0,1),'LineWidth',3);
    %  %hold off;
    %  %colormap(((circshift(lipmanya(16),8,1))));
    %  %caxis(max(abs(caxis))*[-1 1]);

    %  % We assume that T(i,j) is the second vertex of all of these:
    %  G1i = G(J(FMAP(i,mod(j+1-1,4)+1),1),:);
    %  G2i = G(J(FMAP(i,mod(j+2-1,4)+1),1),:);
    %  G3i = G(J(FMAP(i,mod(j+3-1,4)+1),1),:);
    %  M1i = MF(i,mod(j+1-1,4)+1);
    %  M2i = MF(i,mod(j+2-1,4)+1);
    %  M3i = MF(i,mod(j+3-1,4)+1);
    %  G1i(~M1i,:) = fliplr(G1i(~M1i,:));
    %  G2i(~M2i,:) = fliplr(G2i(~M2i,:));
    %  G3i(~M3i,:) = fliplr(G3i(~M3i,:));
    %  assert(G1i(2) == T(i,j) && G2i(2) == T(i,j) && G3i(2) == T(i,j));
    %  assert(G1i(1) ~= G2i(1) && G2i(1) ~= G3i(1) && G3i(1) ~= G1i(1));
    %  TGi = [T(i,j) G1i(1) G2i(1) G3i(1)];
    %  % flip based on j to make volume match sign.
    %  TGi(mod(j,2)==0,:) = TGi(mod(j,2)==0,[2 3 4 1]);
    %  %fprintf('%d %d %d %d %d\n',[M1i M2i M3i sign(volume(V,T(i,:))) == sign(volume(U,TGi)) j]);
    %  %B(end+1,:) = [M1i M2i M3i sign(volume(V,T(i,:))) == sign(volume(U,TGi)) j];
    %  assert(sign(volume(V,T(i,:))) == sign(volume(U,TGi)));
    %  % claim: can determine which orientation to use based on MF(i,:)
    %  TG3 = [TG3;TGi];

    %  %drawnow;
    %end
    %TG3

    i = find(T3);
    [~,j] = min(TF_was_split(i,:),[],2);
    G1 = G(J(FMAP(sub2ind(size(FMAP),i,mod(j+1-1,4)+1)),1),:);
    G2 = G(J(FMAP(sub2ind(size(FMAP),i,mod(j+2-1,4)+1)),1),:);
    G3 = G(J(FMAP(sub2ind(size(FMAP),i,mod(j+3-1,4)+1)),1),:);
    M1 = MF(sub2ind(size(MF),i,mod(j+1-1,4)+1));
    M2 = MF(sub2ind(size(MF),i,mod(j+2-1,4)+1));
    M3 = MF(sub2ind(size(MF),i,mod(j+3-1,4)+1));
    G1(~M1,:) = fliplr(G1(~M1,:));
    G2(~M2,:) = fliplr(G2(~M2,:));
    G3(~M3,:) = fliplr(G3(~M3,:));
    assert(all([G1(:,2) G2(:,2) G3(:,2)]==T(sub2ind(size(T),i,j)),'all'));
    TG3 = [T(sub2ind(size(T),i,j)) G1(:,1) G2(:,1) G3(:,1)];
    TG3(mod(j,2)==0,:) = TG3(mod(j,2)==0,[2 3 4 1]);
    assert(all(sign(volume(V,T(i,:))) == sign(volume(U,TG3))));


    %sB = sortrows(unique(B,'rows'))
    %size(sB)

    codes = [
      0 0 1
      0 1 0
      0 1 1
      1 0 0
      1 0 1
      1 1 0];
    codes = [codes zeros(size(codes,1),1);codes ones(size(codes,1),1)];
    %triplets  = cell(1,1,size(codes,1));
    triplets = cat(3, ...
      [1 2 4 5; 1 2 5 3; 1 4 6 5], ...
      [1 2 4 3; 1 3 4 6; 3 4 6 5], ...
      [1 2 4 3; 1 3 4 5; 1 4 6 5], ...
      [1 2 6 3; 2 3 5 6; 2 4 6 5], ...
      [1 2 5 3; 1 2 6 5; 2 4 6 5], ...
      [1 2 6 3; 2 3 4 6; 3 4 6 5], ...
      [1 2 3 4; 1 2 4 5; 2 4 5 6], ...
      [1 2 3 6; 1 3 4 6; 1 4 5 6], ...
      [1 2 3 4; 1 2 4 6; 1 4 5 6], ...
      [1 2 3 5; 2 3 5 6; 3 4 5 6], ...
      [1 2 3 5; 2 3 5 4; 2 4 5 6], ...
      [1 2 3 6; 1 3 5 6; 3 4 5 6]);
    %TQ3 = [];
    %for i = 1:size(T,1)
    %  if ~T3(i)
    %    continue;
    %  end
    %  % index of vertex on the "cap"
    %  [~,j] = min(TF_was_split(i,:),[],2);
    %  %% These triangles always share an edge with face opposite the cap vertex
    %  %A1i = G(J(FMAP(i,mod(j+1-1,4)+1),2),:);
    %  %A2i = G(J(FMAP(i,mod(j+2-1,4)+1),2),:);
    %  %A3i = G(J(FMAP(i,mod(j+3-1,4)+1),2),:);
    %  %% If you view the tet's face from the exterior:
    %  %%   o            o
    %  %%   |\          /|
    %  %%   |_\        /_|
    %  %%   |B/\      /\B| 
    %  %%   |/ A\    / A\|
    %  %%   o----o  o----o
    %  %%     M=1     M=0
    %  %%
    %  M1i = MF(sub2ind(size(MF),i,mod(j+1-1,4)+1));
    %  M2i = MF(sub2ind(size(MF),i,mod(j+2-1,4)+1));
    %  M3i = MF(sub2ind(size(MF),i,mod(j+3-1,4)+1));
    %  %[[A1i;A2i;A3i] [M1i M2i M3i]'];
    %  %base2 = [A1i(M1i*1+2) A2i(M2i*1+2) A3i(M3i*1+2)];
    %  base = [ ...
    %    G(J(FMAP(i,mod(j+1-1,4)+1),2),M1i*1+2) ...
    %    G(J(FMAP(i,mod(j+2-1,4)+1),2),M2i*1+2) ...
    %    G(J(FMAP(i,mod(j+3-1,4)+1),2),M3i*1+2)];
    %  %assert(all(base == base2));

    %  %% These triangles always share an edge with the new "cap tet"
    %  %B1i = G(J(FMAP(i,mod(j+1-1,4)+1),3),:);
    %  %B2i = G(J(FMAP(i,mod(j+2-1,4)+1),3),:);
    %  %B3i = G(J(FMAP(i,mod(j+3-1,4)+1),3),:);
    %  %%[[B1i;B2i;B3i] [M1i M2i M3i]']
    %  %top2 = [B1i(M1i*1+2) B2i(M2i*1+2) B3i(M3i*1+2)];
    %  top = [ ...
    %    G(J(FMAP(i,mod(j+1-1,4)+1),3),M1i*1+2) ...
    %    G(J(FMAP(i,mod(j+2-1,4)+1),3),M2i*1+2) ...
    %    G(J(FMAP(i,mod(j+3-1,4)+1),3),M3i*1+2)];
    %  %assert(all(top == top2));
    %  code = [M1i M2i M3i mod(j,2)];
    %  [~,code_ind] = ismember(code,codes,'rows');

    %  Ii = [top base];
    %  %Mi = [M1i M2i M3i M1i M2i M3i mod(j+1,2) mod(j,2)];
    %  %Fi = [A1i;A2i;A3i;B1i;B2i;B3i;top;base];
    %  %Fi(Mi==0,:) = fliplr(Fi(Mi==0,:));
    %  %[~,vec_vol] = centroid(U,Fi);
    %  %assert(sign(volume(V,T(i,:))) == sign(vec_vol));

    %  %% reindex Fi to use indices into Ii
    %  %[~,Fi] = ismember(Fi,Ii);
    %  %[~,Ti] = tetrahedralize(U(Ii,:),Fi,'Flags','YY');
    %  %assert(abs(sum(volume(U(Ii,:),Ti))-sum(volume(U(Ii,:),canonicalize_tets(Ti))))<1e-8);
    %  %Ti = canonicalize_tets(Ti);
    %  %cTi = triplets(:,:,code_ind);
    %  %assert(all(cTi == Ti,'all'));
    %  %Ti = canonicalize_tets(Ti);
    %  %if isempty(triplets{code_ind}) || ...
    %  %    any(triplets{code_ind} ~= Ti,'all')
    %  %  triplets{code_ind} = cat(3,triplets{code_ind},Ti);
    %  %end
    %  %Ti = triplets(:,:,code_ind);
    %  Ti = triplets(:,:,code_ind);
    %  TQ3 = [TQ3;Ii(Ti)];

    %  %clf;
    %  %%tsurf(F,V,'CData',D,falpha(0.5,1),fphong);
    %  %hold on;
    %  %%tsurf(edges(G),U,'VertexIndices',1,'CData',S*D,falpha(0.1,1),fphong);
    %  %tsurf(edges([A1i;A2i;A3i;B1i;B2i;B3i]),U,'VertexIndices',1,falpha(0,0.5));
    %  %tsurf([A1i;A2i;A3i],U,'FaceColor',blue,falpha(0.2,0));
    %  %tsurf([B1i;B2i;B3i],U,'FaceColor',orange,falpha(0.2,0));
    %  %%qvr(barycenter(U(Ii,:),Fi),normalizerow(normals(U(Ii,:),Fi)),'k','LineWidth',2);
    %  %%tsurf(Fi,U(Ii,:),'FaceVertexCData',[0 1 1;0 1 1;0 1 1;1 1 0;1 1 0;1 1 0],falpha(0.0,1),'LineWidth',2,'EdgeColor','r');
    %  %tsurf(edges(Ii(Ti)),U,falpha(0.0,1),'LineWidth',2,'EdgeColor','r');
    %  %%tsurf(edges(cTi),U(Ii,:),falpha(0.0,1),'LineWidth',2,'EdgeColor','g');
    %  %tsurf(boundary_faces(T(i,:)),U,'FaceColor',blue,falpha(0.0,1),'LineWidth',1);
    %  %hold off;
    %  %colormap(((circshift(lipmanya(16),8,1))));
    %  %caxis(max(abs(caxis))*[-1 1]);
    %  %pause
    %end

    % Same as above for TG3
    i = find(T3);
    [~,j] = min(TF_was_split(i,:),[],2);
    M1 = MF(sub2ind(size(MF),i,mod(j+1-1,4)+1));
    M2 = MF(sub2ind(size(MF),i,mod(j+2-1,4)+1));
    M3 = MF(sub2ind(size(MF),i,mod(j+3-1,4)+1));
    base = [ ...
      G(sub2ind(size(G),J(FMAP(sub2ind(size(FMAP),i,mod(j+1-1,4)+1)),2),M1*1+2)) ...
      G(sub2ind(size(G),J(FMAP(sub2ind(size(FMAP),i,mod(j+2-1,4)+1)),2),M2*1+2)) ...
      G(sub2ind(size(G),J(FMAP(sub2ind(size(FMAP),i,mod(j+3-1,4)+1)),2),M3*1+2))];
    top = [ ...
      G(sub2ind(size(G),J(FMAP(sub2ind(size(FMAP),i,mod(j+1-1,4)+1)),3),M1*1+2)) ...
      G(sub2ind(size(G),J(FMAP(sub2ind(size(FMAP),i,mod(j+2-1,4)+1)),3),M2*1+2)) ...
      G(sub2ind(size(G),J(FMAP(sub2ind(size(FMAP),i,mod(j+3-1,4)+1)),3),M3*1+2))];
    code = [M1 M2 M3 mod(j,2)];
    [~,code_ind] = ismember(code,codes,'rows');
    TQ32 = triplets(:,:,code_ind);
    TQ32 = TQ32 + permute(6*(0:size(TQ32,3)-1),[1 3 2]);

    TQ32 = reshape(permute(TQ32,[1 3 2]),numel(code_ind)*3,4);
    % 6 by #T3
    Ii = [top base]';
    TQ32 = Ii(TQ32);
    %assert(all(TQ32 == TQ3,'all'));
    TQ3 = TQ32;

    %assert(all(cellfun(@(x) size(x,3),triplets)==1))
    %save('triplets','triplets');
    TG3 = [TG3;TQ3];


    %  TG4 = [];
    %  bad = 0;
    %  for i = 1:size(T,1)
    %    if ~T4(i)
    %      continue;
    %    end
    %    % index of vertex on the "cap"
    %    j = 1;
    %    M1i = MF(i,1);
    %    M2i = MF(i,2);
    %    M3i = MF(i,3);
    %    M4i = MF(i,4);
    %    % Are we guaranteed that T(i,j) occurs as the second vertex of one of
    %    % these?
    %    % Yes.
    %    C1i = G(J(FMAP(i,1),1),:);
    %    C2i = G(J(FMAP(i,2),1),:);
    %    C3i = G(J(FMAP(i,3),1),:);
    %    C4i = G(J(FMAP(i,4),1),:);
    %    C1i(M1i==0,:) = fliplr(C1i(M1i==0,:));
    %    C2i(M2i==0,:) = fliplr(C2i(M2i==0,:));
    %    C3i(M3i==0,:) = fliplr(C3i(M3i==0,:));
    %    C4i(M4i==0,:) = fliplr(C4i(M4i==0,:));
    %    mid = [C1i([1 3]) C2i([1 3]) C3i([1 3]) C4i([1 3])]';
    %    S = [S;sparse(ones(8,1),mid,1/8,1,size(U,1))*S];
    %    U = [U;mean(U(mid,:),1)];
    %    mid_ind = size(U,1);
    %    %mid = [C1i(1);C2i(3);C3i(1);C4i(3)];
    %    %mid
    %    assert(all(mid > size(V,1)));
    %    A1i = G(J(FMAP(i,1),2),:);
    %    A2i = G(J(FMAP(i,2),2),:);
    %    A3i = G(J(FMAP(i,3),2),:);
    %    A4i = G(J(FMAP(i,4),2),:);
    %    A1i(M1i==0,:) = fliplr(A1i(M1i==0,:));
    %    A2i(M2i==0,:) = fliplr(A2i(M2i==0,:));
    %    A3i(M3i==0,:) = fliplr(A3i(M3i==0,:));
    %    A4i(M4i==0,:) = fliplr(A4i(M4i==0,:));
    %    B1i = G(J(FMAP(i,1),3),:);
    %    B2i = G(J(FMAP(i,2),3),:);
    %    B3i = G(J(FMAP(i,3),3),:);
    %    B4i = G(J(FMAP(i,4),3),:);
    %    B1i(M1i==0,:) = fliplr(B1i(M1i==0,:));
    %    B2i(M2i==0,:) = fliplr(B2i(M2i==0,:));
    %    B3i(M3i==0,:) = fliplr(B3i(M3i==0,:));
    %    B4i(M4i==0,:) = fliplr(B4i(M4i==0,:));
    %    %[G1i(2) G2i(2) G3i(2) G4i(2)]
    %    %[G1i(2) G2i(2) G3i(2) G4i(2)]==T(i,j)
    %    %assert(any([G1i(2) G2i(2) G3i(2) G4i(2)]==T(i,j)))
    %    Gi = [A1i;A2i;A3i;A4i;B1i;B2i;B3i;B4i;C1i;C2i;C3i;C4i];
    %    Ti = [repmat(mid_ind,12,1) Gi];
    %    %[volume(U,T(i,:)) sum(volume(U,Ti))]
    %    TG4 = [TG4;Ti];
    %    
    %    %[Uii,~,~,Gii] = remove_unreferenced(U,Gi,true);
    %    %[TUi,TTi] = tetrahedralize(Uii,Gii,'Flags','YY');
    %    %bad = bad + (size(TUi,1) ~= 8);
    %    %
    %    % I'd like to avoid adding tons of Steiner points.
    %    %
    %    % This is happening upwards of 34% of the time on random meshes.
    %    % Even if we try to process tets one by one, _before_ assigning all face
    %    % diagonals, then I'm still not sure we won't get cycles and still have to
    %    % insert Steiner points.
    %    %
    %    % If we went the one-by-one route, we might want to first know how many 
    %    % T4 tets are sharing faces. Each "connected component" of these can be
    %    % processed (in terms of prescribing face diagonals) independently. This
    %    % is surely something needed to be done outside MATLAB.
    %    % 
    %    % If we did one-by-one, I would gather all the T4 tets and then process
    %    % them in descending order of incidence on other T4 tets. That way, a T4
    %    % tet with 4 T4 neighbors gets full choice of face diagonals, with less
    %    % risk of getting trapped later.
    %    %
    %    % The most obvious Steiner point to insert is the "center" of the
    %    % cut plane quad. Then we can just connect all faces to this new point.
    %    %
    %    % If we wanted to keep things parallel. We could still try to identify the
    %    % non-Steiner cases. For these we consider the bottom and top wedges. They
    %    % each can either impose a diagonal on the cut-plane quad or say they
    %    % don't care. If they both impose a choice and they are different, then
    %    % we must add a Steiner point.
    %    %if 8 ~= size(TUi,1)
    %    %  [size(TUi,1)]
    %      %clf;
    %      %%tsurf(F,V,'CData',D,falpha(0.5,1),fphong);
    %      %hold on;
    %      %%tsurf(Gii,Uii,'VertexIndices',1,'FaceVertexCData',[repmat(blue,4,1);repmat(orange,4,1);repmat(0.5*orange,4,1)],falpha(1.0,1));
    %      %tsurf(Gi,U,'VertexIndices',1,'CData',S*D,falpha(0.2,1),fphong);
    %      %qvr(barycenter(U,Gi),normalizerow(normals(U,Gi)),'k','LineWidth',2);
    %      %%tsurf(edges(T(i,:)),U,'CData',S*D,falpha(0.2,1),'LineWidth',2,fphong);
    %      %%tsurf(boundary_faces(T(i,:)),U,'CData',S*D,falpha(0.2,1),'LineWidth',2,fphong);
    %      %sct(U(mid,:),'filled','k','SizeData',500);
    %      %sct(U(end,:),'filled','r','SizeData',500);
    %      %hold off;
    %      %colormap(((circshift(lipmanya(16),8,1))));
    %      %caxis(max(abs(caxis))*[-1 1]);
    %      %axis equal;
    %      %view(13,35);
    %      %pause
    %    %end
    %  end

    % Similar to above
    i = find(T4);
    M1 = MF(i,1);
    M2 = MF(i,2);
    M3 = MF(i,3);
    M4 = MF(i,4);
    C1 = G(J(FMAP(i,1),1),:);
    C2 = G(J(FMAP(i,2),1),:);
    C3 = G(J(FMAP(i,3),1),:);
    C4 = G(J(FMAP(i,4),1),:);
    C1(M1==0,:) = fliplr(C1(M1==0,:));
    C2(M2==0,:) = fliplr(C2(M2==0,:));
    C3(M3==0,:) = fliplr(C3(M3==0,:));
    C4(M4==0,:) = fliplr(C4(M4==0,:));
    A1 = G(J(FMAP(i,1),2),:);
    A2 = G(J(FMAP(i,2),2),:);
    A3 = G(J(FMAP(i,3),2),:);
    A4 = G(J(FMAP(i,4),2),:);
    A1(M1==0,:) = fliplr(A1(M1==0,:));
    A2(M2==0,:) = fliplr(A2(M2==0,:));
    A3(M3==0,:) = fliplr(A3(M3==0,:));
    A4(M4==0,:) = fliplr(A4(M4==0,:));
    B1 = G(J(FMAP(i,1),3),:);
    B2 = G(J(FMAP(i,2),3),:);
    B3 = G(J(FMAP(i,3),3),:);
    B4 = G(J(FMAP(i,4),3),:);
    B1(M1==0,:) = fliplr(B1(M1==0,:));
    B2(M2==0,:) = fliplr(B2(M2==0,:));
    B3(M3==0,:) = fliplr(B3(M3==0,:));
    B4(M4==0,:) = fliplr(B4(M4==0,:));
    mid = size(U,1)+(1:numel(i))';
    G4 = [A1;A2;A3;A4;B1;B2;B3;B4;C1;C2;C3;C4];
    TG4 = [repmat(mid,12,1) G4];
    S4 = sparse( ...
      repmat(1:numel(i),1,8)', ...
      [C1(:,[1 3]) C2(:,[1 3]) C3(:,[1 3]) C4(:,[1 3])], ...
      1/8, ...
      numel(i), ...
      size(U,1))*S;
    S = [S;S4];
    %mid = cat(3,C1(:,1),C1(:,3),C2(:,1),C2(:,3),C3(:,1),C3(:,3),C4(:,1),C4(:,3));
    U4 = mean(cat(3, ...
      U(C1(:,1),:), ...
      U(C1(:,3),:), ...
      U(C2(:,1),:), ...
      U(C2(:,3),:), ...
      U(C3(:,1),:), ...
      U(C3(:,3),:), ...
      U(C4(:,1),:), ...
      U(C4(:,3),:)),3);
    U = [U;U4];



    %nnz(T4)/(nnz(T3)+nnz(T4))
    %bad/nnz(T4)
    TG = [TG0;TG3;TG4];
    I0 = find(T0);
    I3 = find(T3);
    %I3 = [I3;nan

    I4 = find(T4);
    I3 = [I3;repelem(I3,3)];
    I4 = repmat(I4,12,1);
    I = [I0;I3;I4];

    %statistics(U,TG)
    %tsurf(boundary_faces(TG),U,'CData',S*D,falpha(1.0,1),fphong);
    %hold on;
    %tsurf(LE,U,'LineWidth',2);
    %hold off;
    %axis equal;
    %caxis(max(abs(caxis))*[-1 1]);



    G = TG;

    % cheapskate way to get LE
    % Could cull away those not in T3 or T4
    allF =  [...
      G(:,[4 2 3]); ...
      G(:,[3 1 4]); ...
      G(:,[2 4 1]); ...
      G(:,[1 3 2])];
    nv = size(V,1);
    Q = allF>nv;
    N = sum(Q,2) == 3;
    O = G(:);
    assert(all(O(N)<=nv));
    % Avoid computing S*D
    DN = D(O(N));
    LE = allF(N,:);
    LE = LE(DN<0,:);

  case 3
    assert(all(D~=0,'all'),'Level set must not cross exactly at vertices');
    % Find all edge crossings
    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    [E,~,EMAP] = unique(sort(allE,2),'rows');

    % most robust way to check if crossing?
    K = sign(D(E(:,1))) ~= sign(D(E(:,2)));
    % Find the crossing point
    T = -D(E(K,1))./(D(E(K,2))-D(E(K,1)));
    assert(all(T>0 & T<1),'Crossing point must be on edge');
    n = size(V,1);
    nk = nnz(K);
    J = n+(1:nk);
    U = [V;(1-T).*V(E(K,1),:)+T.*V(E(K,2),:)];

    % For every triangle, determine if it's incident on an edge crossing.
    W = reshape(K(EMAP),size(F));
    assert(all( sum(W,2) == 0 | sum(W,2) == 2),'Each triangle can be incident on at most one crossing');
    hot_mask = any(W,2);
    hot = find(hot_mask);
    % find corner that's not across from a crossing.
    [~,C] = min(W(hot,:),[],2);
    E2J = nan(size(E,1),1);
    E2J(K) = n+(1:nk);

    I1 = F(sub2ind(size(F),hot,mod(C+0-1,3)+1));
    I2 = F(sub2ind(size(F),hot,mod(C+1-1,3)+1));
    I3 = F(sub2ind(size(F),hot,mod(C+2-1,3)+1));

    I31 = E2J(EMAP(sub2ind(size(F),hot,mod(C-1-1,3)+1)));
    I12 = E2J(EMAP(sub2ind(size(F),hot,mod(C+1-1,3)+1)));

    not_hot = find(~hot_mask);
    G = [F(~hot_mask,:);I12 I1 I31;I12 I2 I3;I2 I12 I31];
    I = [not_hot;hot;hot;hot];
    nh = numel(not_hot);
    h = numel(hot);
    H = [repmat(0,nh,1);repmat(1,h,1);repmat(2,h,1);repmat(3,h,1)];
    S = [speye(n);sparse([1:nk;1:nk]',E(K,:),[1-T T],nk,n)];

    LE = [I12 I31];
    backwards = D(I1)<0;
    LE(backwards,:) = LE(backwards,[2 1]);
  end

end
