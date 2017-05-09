function [SV,SF,J,DT] = linear_sweep(V,F,sweep,varargin)
  % LINEAR_SWEEP Compute the surface of the solid sweep of a surface mesh (V,F)
  % along a vector (sweep): i.e. the Minkowski sum of (V,F) and the line
  % segment (0,0)-->(sweep) 
  %
  % [SV,SF] = linear_sweep(V,F,sweep)
  % [SV,SF] = linear_sweep(V,F,sweep,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by dim list of facet indices into V
  %   sweep  1 by dim vector defining sweep direction and distance
  %   Optional:
  %     'SelfUnion' followed by whether to resolve intersections and extract
  %       the self union of the result. Only relevant in 3D. Faster if `false`
  %       but result will be self-intersecting for non-convex inputs {true}.
  % Outputs:
  %   SV  #SV by dim list of vertex positions
  %   SF  #F by dim list of facet indices into SV
  %   J  #F list of indices into [F;F+v;sweep] revealing birth parents 
  %

  % default values
  self_union = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'SelfUnion'}, ...
    {'self_union'});
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
    SE = fliplr(outline(TF(round(w)~=0,:)));
    [SV,IM] = remove_unreferenced(TV,SE);
    SE = IM(SE);
    % output names
    SF =SE;
  case 3
    % Flip so all are facing along with v
    B = F;
    % number of triangles
    m = size(F,1);
    %% unzip mesh
    %W = V(B,:);
    %B = [1:m;m+[1:m;m+(1:m)]]';
    W = V;
    B = F;
    N = normals(W,B);
    DT = sum(bsxfun(@times,N,sweep),2);
    if any(DT==0)
      warning('faces parallel to sweep not well supported');
    end
    A = DT>0;
    BA = B;
    BA(A,:) = fliplr(BA(A,:));
    % number of vertices
    n = size(W,1);
    % duplicate
    W = [W;bsxfun(@plus,W,sweep)];
    %% Prism for each triangle
    %GT = [ ...
    %  % Original source
    %  B; ...
    %  % Original dest
    %  n+B; ...
    %  % Bottom
    %  BA; ...
    %  % Top
    %  n+fliplr(BA); ...
    %  ];
    %[S,GT] = total_signed_occurrences(GT);
    %GT = [GT(S==2,:);fliplr(GT(S==-2,:))];
    % Or equivalently 
    JB = (1:m)';
    GT = [B(~A,:); n+B(A,:);];
    JT = [JB(~A); m+JB(A);];
    % This is more general for the case that (V,F) is already non-manifold
    GQ = [ ...
      % Sides
      BA(:,[2 1]) n+BA(:,[1 2]); ...
      BA(:,[3 2]) n+BA(:,[2 3]); ...
      BA(:,[1 3]) n+BA(:,[3 1]); ...
      ];
    [S,GQ] = total_signed_occurrences(GQ);
    GQ = [GQ(S==2,:);fliplr(GQ(S==-2,:))];
    % triangulate
    G = [GT;GQ(:,[1 2 3]);GQ(:,[1 3 4])];
    J = [JT;2*m+ones(2*size(GQ,1),1)];
    if self_union
      [SV,SF,SJ] = mesh_boolean(W,G,[],[],'union');
      J = J(SJ);
    else
      SV = W;
      SF = G;
    end
  end

end
