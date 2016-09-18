function [U,G,J,BC,SU,L] = slice_isolines(V,F,SV,val,varargin)
  % SLICE_ISOLINES Slice through a triangle mesh (V,F) at isolines val of a given
  % per-vertex scalar function.
  %
  % [U,G] = slice_isolines(V,F,SV,val);
  % [U,G,J,BC,SU] = slice_isolines(V,F,SV,val, ...
  %   'ParameterName',parameter_value, ...)
  %
  % Inputs:
  %   V  #V by dim list of tet mesh vertices
  %   F  #F by 3 list of tet indices into V 
  %   S  #V list of scalar values per vertex
  %   val  #val list of isolines values
  %   Optional:
  %     'Manifold' followed by whether to stitch together triangles into a
  %       manifold mesh {false}: results in more compact U but slightly slower.
  % Outputs:
  %   U  #U by 3 list of triangle mesh vertices along slice
  %   G  #G by 3 list of triangles indices into U
  %   J  #G list of indices into F revealing which tet this face came from
  %   BC  #U by #V list of barycentric coordinates (or more generally: linear
  %     interpolation coordinates) so that U = BC*V
  %   SU  #U list of interpolated scalar values at U
  %   L  #G list of indices into val of next value above this triangle
  %

  manifold = false;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Manifold'}, ...
    {'manifold'});
  v = 1;
  while v <= numel(varargin)
    param_name = varargin{v};
    if isKey(params_to_variables,param_name)
      assert(v+1<=numel(varargin));
      v = v+1;
      % Trick: use feval on anonymous function to use assignin to this
      % workspace
      feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
    else
      error('Unsupported parameter: %s',varargin{v});
    end
    v=v+1;
  end

  % helper assuming val is a scalar
  function [U,G,J,BC,SU,L] = single_val(V,F,SV,val)
    function [U,G,BC,L] = one_below(V,F,SF,val)
      [sSF,sJ] = sort(SF,2);
      sF = F(sub2ind(size(F),repmat(1:size(F,1),size(F,2),1)',sJ));
      lambda = (sSF(:,2:3)-val)./bsxfun(@minus,sSF(:,2:3),sSF(:,1));
      BC = sparse( ...
        repmat((1:size(sF,1)*2)',1,2), ...
        [repmat(sF(:,1),2,1) reshape(sF(:,2:3),size(sF,1)*2,1)], ...
        [lambda(:) 1-lambda(:)], ...
        size(sF,1)*2,size(V,1));
      U = [V;BC * V];
      % Split into three triangles
      G = [ ...
        sF(:,1) size(V,1)+[1:size(F,1);size(F,1)+(1:size(F,1))]'; ...
        fliplr([sF(:,2) size(V,1)+[1:size(F,1);size(F,1)+(1:size(F,1))]']); ...
        sF(:,[2 3]) size(V,1)+[size(F,1)+(1:size(F,1))]'; ...
        ];
      L = [ones(size(sF,1),1);2*ones(2*size(sF,1),1)];
      flip = repmat( ...
        (sJ(:,1)==1 & sJ(:,2)==3) | ...
        (sJ(:,1)==2 & sJ(:,2)==1) | ...
        (sJ(:,1)==3 & sJ(:,2)==2),3,1);
      G(flip,:) = fliplr(G(flip,:));
    end
    SF = SV(F);
    I12 = sum(SF<val,2) == 1;
    % U is the running set of vertices
    U = V;
    [U,G12,BC12,L12] = one_below(U,F(I12,:),SF(I12,:),val);
    I21 = sum(SF>=val,2) == 1;
    [U,G21,BC21,L21] = one_below(U,F(I21,:),2*val-SF(I21,:),val);
    BC = [speye(size(V,1));BC12;BC21(:,1:size(V,1))];
    untouched = find(~I12&~I21);
    G = [F(untouched,:);G12;G21];
    J = [untouched;repmat(find(I12),3,1);repmat(find(I21),3,1)];
    SU = BC*SV;
    L = [(sum(SF(untouched,:)>=val,2)==3)+1;L12;3-L21];
  end

  U = V;
  G = F;
  SU = SV;
  J = 1:size(F,1);
  L = ones(size(F,1),1);
  for v = 1:numel(val)
    prev_J = J;
    prev_L = L;
    [U,G,J,BC,SU,L] = single_val(U,G,SU,val(v));
    L = prev_L(J)+(L-1);
    J = prev_J(J);
  end

  if manifold
    % should be able to do this combinatorially
    bbd = normrow(max(V)-min(V));
    flag = [(1:size(V,1))';ones(size(U,1)-size(V,1),1)];
    [U,I,IM] = remove_duplicate_vertices([U flag],1e-14*bbd);
    U = U(:,1:3);
    BC = BC(I,:);
    SU = SU(I,:);
    G = IM(G);
  end

end
