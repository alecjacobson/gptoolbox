function [VH] = inflate(varargin)
  %INFLATE a planar 2D mesh into a 3D "pillow" with the original silhouette 
  %
  % [VH] = inflate(V,F,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by 2 list of mesh vertices
  %   F  #F by 3 list of triangle indices into
  %  Optional:
  %    'Method' followed by one of the following: {'BiXYThenBaran'}
  %      'Baran' Follow "Notes on Inflating Curves" [Baran and Lehtinen 2009].
  %        Leave X and Y coordinates as is and solve Δ(z^2) = f, that is solve
  %        Δw = f and set z = sqrt(w), 
  %      'BiXYHarmonic', solve Δ²XY = 0 with implicit neumann conditions and
  %        simultaneous solve Δz = f, 
  %      'BiXYThenBaran', solve Δ²XY = 0 with implicit neumann conditions,
  %        stop, recompute cotmatrix and solve z as in 'Baran', 
  %      'MixedFEM' solve Δ²XYZ = 0 with ∂XY/∂n = 0 and ∂Z/∂n = f, 
  %      'Repousse' follow "Repousse: Automatic Inflation of 2D Artwork"
  %        [Joshi and Carr 2008], that is solve Δ²XY = 0 with ∂XYZ/∂n = 0, but
  %        explicitly treat ΔZ|∂Ω = f/2
  %    'PuffFactor' followed by a scale value corresponding to different amount of
  %      "puffiness". See f value in 'Method' descriptions {f=-8}
  % Outputs:
  %   VH  #V by 3 list of new "height" field vertices, but note that depending
  %     on 'Method', XY values may be different than V
  %
  % See also: triangle, outline, get_pencil_curve, glue_reverse
  % 

  V = varargin{1};
  F = varargin{2};
  % get outline ("boundary") of mesh 
  out = unique(reshape(outline(F),[],1));
  no = numel(out);

  % default optional parameter values
  method = 'BiXYThenBaran';
  f = -8;

  ii = 3;
  while(ii <= nargin)
    switch varargin{ii}
    case 'Method'
      ii = ii + 1;
      assert(ii<=nargin);
      method = varargin{ii};
    case 'PuffFactor'
      ii = ii + 1;
      assert(ii<=nargin);
      f = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii + 1;
  end

  % mesh-independence scale
  scale = 2/max((max(V)-min(V)));
  % number of domain vertices
  n = size(V,1);
  L = cotmatrix(V,F);
  M = massmatrix(V,F,'voronoi');

  switch method
  case 'Baran'
    Z2 = min_quad_with_fixed(-L,f*sqrt(scale)*M*ones(n,1),out,zeros(no,1));
    VH = [V sqrt(Z2)];
  case 'BiXYHarmonic'
    % solve poisson equation for Z-coord
    Z = min_quad_with_fixed(-L,f*scale*M*ones(n,1),out,zeros(no,1));
    % solve bipoisson equation for XY-coords
    VH = min_quad_with_fixed(L*(M\L),zeros(n,1),out,V(out,1:2));
    VH = [VH Z];
  case 'BiXYThenBaran'
    % solve bipoisson equation for XY-coords
    VH = min_quad_with_fixed(L*(M\L),zeros(n,1),out,V(out,1:2));
    % recompute cotmatrix using VH
    L = cotmatrix(VH,F);
    M = massmatrix(VH,F,'voronoi');
    %tsurf(F,VH)
    %input('')
    Z2 = min_quad_with_fixed(-L,f*sqrt(scale)*M*ones(n,1),out,zeros(no,1));
    VH = [VH sqrt(Z2)];
  case 'MixedFEM'
    % "Mixed Finite Elements for Variational Surface Modeling" [Jacobson et al.
    % 2010]
    indices = 1:n;
    Omega = indices(~ismember(indices,out));
    N0 = out(:)';
    all = [N0 Omega];
    bn = zeros(n,3);
    bn(out,3) = 1/scale/-f;
    bn = bn(all,:);
    b0 = [V(N0,1:2) zeros(numel(N0),1)];
    flatten = false;
    if flatten
      A = L(Omega,all) * (M(all,all) \ L(all,Omega));
      %rhs_Dx = -L(all,N0)*b0 +  bn;
      %rhs = L(Omega,all) * (M(all,all) \ rhs_Dx);
      %rhs = L(Omega,all) * (M(all,all) \ ...
      %  (-L(all,N0)*b0)) +  L(Omega,all) * (M(all,all) \ bn);
      rhs = ...
        L(Omega,all) * (M(all,all) \ (-L(all,N0)*b0)) +  ...
        L(Omega,all) * (M(all,all) \ bn);
    else
      Z = sparse(n,n);
      A = [-M(  all,all) L(  all,Omega);
            L(Omega,all) Z(Omega,Omega)];
      rhs_Dx = -L(all,N0)*b0 + bn;
      rhs_Dy = zeros(numel(Omega),size(rhs_Dx,2));
      rhs = [rhs_Dx;rhs_Dy];
    end
    sol = A\rhs;
    VH = [V(:,1:2) zeros(n,1)];
    VH(Omega,:) = sol((end-numel(Omega)+1):end,:);
  case 'Repousse'
    % "Repousse??: Automatic Inflation of 2D Artwork" [Joshi and Carr 2008]
    Z = sparse(n,n);
    % dirichlet on boundary and fixed laplacian on boundary
    b = [out;n+out];
    bc = [[V(out,1:2) zeros(no,1)]; repmat([0,0,f/2*scale],no,1)];
    VH = min_quad_with_fixed([Z -L;-L M],zeros(2*n,1),b,bc);
    % only keep positions 
    VH = VH(1:n,:);
  otherwise
    error(['Unsupported method: ' varargin{ii}]);
  end


end
