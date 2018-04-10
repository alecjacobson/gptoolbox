function [CV,CE,I,CN,GD,D] = polygonize(V,F,fun,varargin)
  % Polygonize (contour) an implicit function in the spirit of "An Implicit
  % Surface Polygonizer" [Bloomenthal 1994]
  %
  % [CV,CE] = polygonize(V,F,fun)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions of original mesh
  %   F  #F by dim+1 list of element indices into V
  %   fun  function handle so that zero is the desired level set
  % Outputs:
  %   CV  #CV by dim list of contour mesh vertices 
  %   CE  #CE by dim list of facet indices into CV
  %   I  #CE list of indices into F
  %   CN  #CN by dim list of unit normal vectors at vertices
  %

  % default values
  delta = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Delta'}, ...
    {'delta'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous funtion to use assignin to this workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  flipped_order = flipped_tet_orders();

  max_iters = 10;
  D = fun(V);
  interval = @(DF) any(DF>0.0,2) & any(DF<=0.0,2);
  crosses = @(DE) ...
    (DE(:,1)>0 & DE(:,2)<= 0)*-1 + ...
    (DE(:,2)>0 & DE(:,1)<= 0)*1;
  crossing_F = find(interval(D(F)));
  FF = F(crossing_F,:);
  % simplex size
  ss = size(F,2);
  switch ss
  case 3
    allE = [FF(:,[2 3]);FF(:,[3 1]);FF(:,[1 2])];
  case 4
    allE = ...
      [FF(:,1) FF(:,2); ...
       FF(:,1) FF(:,3); ...
       FF(:,1) FF(:,4); ...
       FF(:,2) FF(:,3); ...
       FF(:,2) FF(:,4); ...
       FF(:,3) FF(:,4) ...
       ];
  end
  sE = sort(allE,2);
  [E,~,EMAP] = unique(sE,'rows');
  cross_dir = crosses(D(E));
  crossing = cross_dir ~= 0;
  J = (1:size(E,1))';
  EE = E(crossing,:);
  J(crossing) = 1:size(EE,1);
  % Blasphemy
  switch ss
  case 3
    CE = sort(reshape(crossing(EMAP),[],3).*reshape(J(EMAP),[],3),2);
    CE = CE(:,2:end);
  case 4
    % CE(f,i) = 0 if ith edge of element f does not cross, otherwise
    % CE(f,i) is the index of the unique edge that does cross
    CE = reshape(crossing(EMAP),[],6).*reshape(J(EMAP),[],6);
    % If 3 edges cross then we can surface with a single triangle
    crossE = reshape(cross_dir(EMAP),[],6);
    IT = find(sum(CE>0,2)==3);
    crossT = crossE(IT,:);
    [CT,ST] = sort(CE(IT,:),2);
    CT = CT(:,end-2:end);
    flippedT = [
      -1 -1 -1 0 0 0 4 5 6 1 3 2
      -1 -1 -1 0 0 0 4 5 6 2 1 3
      -1 0 0 1 1 0 2 3 6 1 5 4
      0 -1 0 -1 0 -1 1 3 5 2 4 6
      0 -1 0 -1 0 1 1 3 5 2 4 6
      0 -1 0 1 0 1 1 3 5 2 4 6
      0 0 -1 0 -1 -1 1 2 4 3 6 5
      0 0 -1 0 1 -1 1 2 4 3 6 5
      0 0 1 0 1 -1 1 2 4 3 5 6
      0 0 1 0 1 1 1 2 4 3 5 6
      0 1 0 -1 0 -1 1 3 5 2 6 4
      1 0 0 -1 -1 0 2 3 6 1 4 5
      1 0 0 1 -1 0 2 3 6 1 4 5
      1 0 0 1 1 0 2 3 6 1 4 5
      1 1 1 0 0 0 4 5 6 1 2 3
      1 1 1 0 0 0 4 5 6 2 3 1
    ];
    fT = ismember([crossT ST],flippedT,'rows');
    CT(fT,:) = fliplr(CT(fT,:));

    IQ = find(sum(CE>0,2)==4);
    crossQ = crossE(IQ,:);
    [CQ,SQ] = sort(CE(IQ,:),2);
    CQ = CQ(:,end-3:end);
    flippedQ = [
      -1 -1 0 0 -1 1 3 4 2 1 6 5
      -1 -1 0 0 1 1 3 4 2 1 6 5
      -1 0 -1 -1 0 -1 2 5 1 3 4 6
      -1 0 -1 1 0 -1 2 5 1 3 4 6
      -1 0 -1 1 0 1 2 5 1 3 4 6
      0 -1 -1 -1 -1 0 1 6 3 2 5 4
      0 1 1 -1 -1 0 1 6 2 3 4 5
      0 1 1 -1 1 0 1 6 2 3 4 5
      0 1 1 1 1 0 1 6 2 3 4 5
      1 0 1 1 0 1 2 5 3 1 6 4
      1 1 0 0 -1 -1 3 4 1 2 5 6
      1 1 0 0 -1 1 3 4 1 2 5 6
    ]; 
    fQ = ismember([crossQ SQ],flippedQ,'rows');
    %% huh. Need to flip _after_ creating triangles.
    %CQT([fQ;fQ],:) = fliplr(CQT([fQ;fQ],:));
    % or switch the middle two...
    CQ(fQ,:) = CQ(fQ,[1 3 2 4]);
    I = [IT;IQ;IQ];
    I = crossing_F(I);
  end
  % Upper and lower bound on barycenteric coordinate locating =0.5
  EEl= zeros(size(EE,1),1);
  EElV = V(EE(:,1),:);
  Dl = D(EE(:,1));
  EEu= ones(size(EE,1),1);
  Du = D(EE(:,2));
  EEuV = V(EE(:,2),:);
  for iter = 1:max_iters
    EEm = (EEl+EEu)/2;
    CV = 0.5*(EElV+EEuV);
    if iter < max_iters
      Dm = fun(CV);
      front = interval([Dl Dm]);
       EEu(front,:) =  EEm(front,:);
      EEuV(front,:) = CV(front,:);
       Du(front,:) =  Dm(front,:);
       EEl(~front,:) =  EEm(~front,:);
      EElV(~front,:) = CV(~front,:);
       Dl(~front,:) =  Dm(~front,:);
    end
  end

  switch ss
  case 4
    % Flip non-deluanay edges in quads... at least.
    CQT = [CQ(:,[3 1 4]);CQ(:,[2 4 1])];
    A = reshape(internalangles(CV,CQT),[],6);
    % Quads with Non-delaunay diagonals.
    nd = (A(:,1)+A(:,4))>pi;
    CQ(nd,:) = CQ(nd,[3 4 1 2]);
    CQ1 = CQ(:,[3 1 4]);
    CQ2 = CQ(:,[2 4 1]);
    CQ1(nd,:) = fliplr(CQ1(nd,:));
    CQ2(nd,:) = fliplr(CQ2(nd,:));
    CQT = [CQ1;CQ2];
    %% Double check: sum should be zero
    % A = reshape(internalangles(CV,CQT),[],6);
    % % Non-delaunay
    % nd = (A(:,1)+A(:,4))>pi;
    % sum(nd)
    CE = fliplr([CT;CQT]);
    assert(size(CE,2) == ss-1);

    %% How I determined what should be flipped:
    % S = sum(normals(CV,CT).*barycenter(CV,CT),2);
    % A = unique([crossT(S<0,:) ST(S<0,:)],'rows');
    % B = unique([crossT(S>0,:) ST(S>0,:)],'rows');
    % intersect(A,B,'rows');
    % fprintf('    flippedT = [\n');
    % fprintf('      %d %d %d %d %d %d %d %d %d %d %d %d\n',A');
    % fprintf('    ];\n');
    %S = sum(normals(CV,CQ(:,[1 4 3])).*barycenter(CV,CQ(:,[1 4 3])),2);
    %A = unique([crossQ(S<0,:) SQ(S<0,:)],'rows');
    %B = unique([crossQ(S>0,:) SQ(S>0,:)],'rows');
    %intersect(A,B,'rows')
    %fprintf('    flippedQ = [\n');
    %fprintf('      %d %d %d %d %d %d %d %d %d %d %d %d\n',A');
    %fprintf('    ];\n');
  end

  if nargout>=4
    if isempty(delta)
      %  following Bloomenthal's "An Implicit Surface Polygonizer" 1994
      %
      delta = avgedge(V,F)/max_iters^2;
      % More accurate...
      delta = delta*0.0001;
    end
    % Approximate the gradient with finite differences (it'd be cool if fun
    % could also provide the gradient)
    %
    % Combine FD calls into a single call to fun in case fun has a lot of
    % per-call precomputation.
    order = 2
    switch order
    case 1
      % Wasn't computed on the last iteration
      D = fun(CV);
      switch ss
      case 3
        GD = ...
          reshape(fun([CV+[1 0]*delta;CV+[0 1]*delta]),[],2);
      case 4
        GD = reshape( ...
          fun([CV+[1 0 0]*delta;CV+[0 1 0]*delta;CV+[0 0 1]*delta]), ...
          size(CV,1),3);
      end
      % Don't bother dividing by delta (we're normalizing anyway)
      G = GD-D;
    case 2
      delta = delta/2;
      switch ss
      case 3
        GD = ...
          reshape(fun([CV+[1 0]*delta;CV+[0 1]*delta;CV-[1 0]*delta;CV-[0 1]*delta]),[],2,2);
      case 4
        GD = reshape( ...
          fun([ ...
            CV+[1 0 0]*delta;CV+[0 1 0]*delta;CV+[0 0 1]*delta; ...
            CV-[1 0 0]*delta;CV-[0 1 0]*delta;CV-[0 0 1]*delta; ...
            ]), ...
          size(CV,1),3,2);
      end
      G = GD(:,:,1)-GD(:,:,2);
    end

    CN = -normalizerow(G);
  end

end
