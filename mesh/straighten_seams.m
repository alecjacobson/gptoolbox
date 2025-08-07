function [UE,UT,OT,SV,OF,C] = straighten_seams(V,F,VT,FT,varargin)
  % STRAIGHTEN_SEAMS Given a obj-style mesh with (V,F) defining the geometric
  % surface of the mesh and (VT,FT) defining the
  % parameterization/texture-mapping of the mesh in the uv-domain, find all
  % seams and boundaries in the texture-mapping and "straighten" them,
  % remapping vertices along the boundary and in the interior. This will be
  % careful to consistently straighten multiple seams in the texture-mesh
  % corresponding to the same edge chains in the surface-mesh. 
  %
  % [UE,UT] = straighten_seams(V,F,VT,FT)
  %
  % Inputs:
  %  V  #V by 3 list of vertices
  %  F  #F by 3 list of triangle indices
  %  VT  #VT by 2 list of texture coordinates
  %  FT  #F by 3 list of triangle texture coordinates
  %  Optional:
  %    'Tol'  followed by Ramer-Douglas-Peucker tolerance as a fraction of the
  %      curves bounding box diagonal (see dpsimplify)
  % Outputs:
  %   UE  #UE by 2 list of indices into VT of coarse output polygon edges
  %   UT  #VT by 3 list of new texture coordinates
  %
  % See also: simplify_curve, dpsimplify

  tol = inf;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Tol'}, ...
    {'tol'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  m = size(F,1);
  assert(m == size(FT,1));
  % boundary edges of the texture-mapping
  [~,BT] = on_boundary(FT);
  [~,BF] = on_boundary(F);
  assert(~any(BF(:) & ~BT(:)), ...
    'Not dealing with boundaries in 3D that get "stitched" during UV mapping');
  ET = [FT(:,[2 3]);FT(:,[3 1]);FT(:,[1 2])];
  % Add "texture occluding contours" and "texture silhouettes", these are
  % combinatorially non-boundary edges but edges with (at least) one positive
  % area triangle and (at least) one negative area triangle (assuming
  % consistent vertex ordering.
  [uET,~,ETMAP] = unique(sort([FT(:,[2:3]);FT(:,[3 1]);FT(:,[1:2])],2),'rows');
  %A = sparse(ETMAP,(1:size(FT,1)*3)',2*(repmat(doublearea(VT,FT),3,1)>0)-1);
  %SIL = min(A,[],2) == -1 & max(A,[],2) == 1;
  %BT = BT | reshape(SIL(ETMAP),size(BT));

  OT = ET(BT(:),:);
  % "half"-edges with indices into 3D-mesh
  EF = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
  OF = EF(BT(:),:);
  % Find unique (undirected) edges in F
  sEF = sort(EF,2);
  if min(sEF(:))>0 
    [~,~,EFMAP] = unique(sEF(:,1)+(max(sEF(:,1)))*sEF(:,2));
  else
    [~,~,EFMAP] = unique(sEF,'rows');
  end
  OFMAP = EFMAP(BT(:));
  % Two boundary edges on the texture-mapping are "equivalent" to each other on
  % the 3D-mesh if their 3D-mesh vertex indices match
  OEQ = sparse((1:size(OT,1))',OFMAP,1,size(OT,1),m*3);
  OEQ = (OEQ*OEQ')~=0;
  OEQ = OEQ-diag(diag(OEQ));
  % For each edge in OT, for each endpoint, how many _other_ texture-vertices
  % are images of all the 3d-mesh vertices in F who map from "corners" in F/FT
  % mapping to this endpoint.
  %
  % Adjacency matrix between 3d-vertices and texture-vertices
  V2VT = sparse(F(:),FT(:),1,size(V,1),size(VT,1))~=0;
  % For each 3d-vertex count how many different texture-coordinates its getting
  % from different incident corners
  DV = full(sum(V2VT,2));
  [M,I] = max(V2VT,[],1);
  assert(all(M==1));
  % Map counts onto texture-vertices
  DT = DV(I);
  % Boundary in 3D && UV
  BTF = BF(BT(:));
  SV = false(size(VT,1),1);
  % Texture-vertex is "sharp" if incident on "half-"edge that is not a boundary
  % in the 3D mesh but is a boundary in the texture-mesh AND is not "cut
  % cleanly" (the vertex is mapped to exactly 2 locations)
  SV(OT(~BTF,:)) = true;
  cuts = sum(sparse(OT,repmat(1:size(OT),2,1)',1,size(VT,1),size(OT,1))*OEQ,2);
  %tsurf(FT,VT,'Cdata',full(cuts==4)*1,fphong);
  %colorbar
  %pause
  CL = DT==2 & cuts==2;
  %CL = CL | (ismember((1:size(VT,1))',uET(SIL,:)) && cuts==4);
  SV(CL) = false;

  % vertices at the corner of ears are declared to be sharp. This is
  % conservative: for example, if the ear is strictly convex and stays strictly
  % convex then the ear won't be flipped.
  [ears,ear_opp] = find_ears(FT);
  earT = sparse(FT(sub2ind(size(FT),ears,ear_opp)),1,1,size(VT,1),1);
  % Even if ear-vertices are marked as sharp if it changes, e.g., from convex
  % to concave then it will _force_ a flip of the ear triangle. So, declare
  % that neighbors of ears are also sharp.
  earT = earT | adjacency_matrix(FT)*earT;

  SV(earT) = true;
  % There might be an sharp vertex on one copy, so mark vertices on other
  % copies, as they live on the 3D mesh
  SV = any(V2VT(:,SV),2);
  SV = any(V2VT(SV,:),1)';
  SV(~ismember(1:end,OT)) = false;

  OTVT = sparse(repmat(1:size(OT,1),2,1)',OT,1,size(OT,1),size(VT,1));
  % Seam components stop at sharp vertices
  A = OTVT*diag(sparse(~SV))*OTVT';
  [nc,C] = conncomp(A);

  %clf;
  %%CC = C;
  %%%CC(~ismember(C,[13 130    95    96    34   158   158   123    95   158])) = nan;
  %for p = 1:2
  %  %subplot(1,2,p);
  %  hold on;
  %  %switch p
  %  %case 2
  %    BC = barycenter(VT,OT);
  %    [OEQI,OEQJ] = find(tril(OEQ));
  %    %BC(isnan(CC),:) = nan;
  %    A = (V2VT'*1*V2VT~=0);
  %    A(:,~SV) = 0;
  %    [AI,AJ] = find(triu(A,1));
  %    %plot_edges(VT,[AI AJ],'-k');
  %    %plot_edges(BC,[OEQI OEQJ],':k');
  %    scatter(VT(SV,1),VT(SV,2),'o','SizeData',300);
  %    tsurf(FT,VT,'FaceAlpha',0.05,'EdgeAlpha',0.02,'FaceColor',blue);
  %    tsurf(reshape([1:size(OT,1),1:size(OT,1)*2],[],3),VT(OT(:),:),'CData',repmat(C',2,1),'EdgeColor','interp','LineWidth',3);
  %    %tsurf(reshape([1:size(OT,1),1:size(OT,1)*2],[],3),VT(OT(:),:),'CData',repmat(CC',2,1),'EdgeColor','interp','LineWidth',3);
  %    %text(BC(:,1),BC(:,2),num2str(C'),'BackgroundColor',[0.5 0.5 0.5]);
  %  %case 1
  %  %  tsurf(F,V,'FaceAlpha',0.8,'EdgeAlpha',0.4);
  %  %  tsurf(reshape([1:size(OF,1),1:size(OF,1)*2],[],3),V(OF(:),:),'CData',repmat(CC',2,1),'EdgeColor','interp','LineWidth',3);
  %  %end
  %  axis equal;
  %  hold off;
  %end
  %error


  % New texture-vertex locations
  UT = VT;
  % Indices into UT of coarse output polygon edges
  UE = [];
  % loop over each component
  done = false(nc,1);
  for c = 1:nc
    % may have handled the copy already
    if done(c)
      continue
    end
    done(c) = true;
    % edges of this component
    Ic = find(C==c);
    if isempty(Ic)
      continue;
    end
    % number of compies
    ncopies = sum(OEQ(Ic(1),:),2)+1;
    % All edges should agree on the number of copies
    assert(all(sum(OEQ(Ic,:),2) == ncopies-1));
    % number of copies might be >2 for non-manifold meshes
    assert(ncopies == 1 || ncopies == 2,'Not dealing with non-manifold meshes');
    switch ncopies
    case 1
      [vpath,epath] = edges_to_path(OT(Ic,:));
      assert(vpath(1) ~= vpath(end) || ~any(SV(vpath)), ...
        'Not dealing with 1-loops touching "sharp" corners');
      % simple open boundary
      PI = VT(vpath,:);
      bbd = norm(max(PI)-min(PI));
      no_boundary_collapse = true;
      assert(numel(vpath) >= 2);
      is_closed = vpath(1) == vpath(end);
      assert(~is_closed || numel(vpath)>=4);
      eff_tol = min(tol,2);
      while true
        [UPI,UIc,UT(vpath,:)] = simplify_curve(PI,eff_tol*bbd);
        if ~is_closed || ~no_boundary_collapse
          break;
        end
        if size(UPI,1)>=4
          break;
        else
          eff_tol = eff_tol/2;
        end
      end
      UE = [UE;reshape(vpath(UIc([1:end-1;2:end]')),[],2)];
    case 2
      % Find copies
      [Icc,II] = find(OEQ(:,Ic));
      assert(isequal(II,(1:numel(Ic))'));
      assert(numel(Icc) == numel(Ic));
      cc = C(Icc(1));
      assert(all(C(Icc) == cc));
      assert(~done(cc));
      done(cc) = true;
      assert(isempty(setxor(OF(Ic,:),OF(Icc,:))));
      % Is copy going the same direction?
      flipped = OF(Ic,1) ~= OF(Icc,1);
      if numel(Ic) == 1
        % No change to UT
        UE = [UE;OT(Ic,:)];
        assert(numel(Icc) == 1);
        if flipped
          UE = [UE;fliplr(OT(Icc,:))];
        else
          UE = [UE;OT(Icc,:)];
        end
      else
        [vpath,epath,eend] = edges_to_path(OT(Ic,:));
        % Flip endpoints if needed
        eend(flipped) = 3-eend(flipped);
        vpathc = ...
          OT(sub2ind(size(OT),Icc(epath([1:end end])),[eend;3-eend(end)]));
        PI = [VT(vpath,:) VT(vpathc,:)];
        bbd = norm(max(PI)-min(PI));
        [~,UIc,SI] = simplify_curve(PI,tol*bbd);
        UT(vpath,:) = SI(:,1:2);
        UT(vpathc,:) = SI(:,3:4);
        UE = [ ...
          UE; ...
          reshape(vpath (UIc([1:end-1;2:end]')),[],2); ...
          reshape(vpathc(UIc([1:end-1;2:end]')),[],2)];
      end
    end
  end
  if nargout > 1
    b = unique(OT);
    %UT = kharmonic(VT,FT,b,UT(b,:),1);
  end
end
