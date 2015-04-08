function [H,L,BT,CC] = arap_hessian(V,F,varargin)
  % ARAP_HESSIAN Compute the Hessian of an as-rigid-as-possible energy defined
  % for a mesh (V,F), according to "Shape Decomposition using Modal Analysis"
  % [Huang et al. 2009].
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by simplex-size list of simplex indices into V
  %   Optional:
  %     'Energy' followed by 'spokes','spokes-and-rims',or {'elements'}
  %     'EvaluationPoint' followed by U  #V by dim list of current positions
  %     'Rotations' followed by best fit rotations R  dim by dim #R
  %       corresponding to U {compute using `fit_rotations`}
  % Outputs:
  %   H  #V*dim by #V*dim 
  %
  % See also: arap, arap_gradient
  % 
  % Example:
  % % Given a mesh (V,F) and deformed positions U0, flow to energy minimum
  % % using Newton's method.
  % clf;
  % hold on;
  %   tsurf(F,V,'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','r');
  %   t = tsurf(F,U0,'FaceAlpha',0.2,'EdgeAlpha',0.2,'FaceColor','b');
  % hold off;
  % axis equal;
  % U = U0;
  % delta_t = 1;
  % while true
  %   [G,E,R] = arap_gradient(V,F,U);
  %   H = arap_hessian(V,F,'EvaluationPoint',U,'Rotations',R);
  %   U = U - reshape(delta_t * (H \ G(:)),size(U));
  %   U = bsxfun(@plus,U,mean(V)-mean(U));
  %   t.Vertices = U;
  %   title(sprintf('E = %g\n',E),'FontSize',20);
  %   drawnow;
  % end
  %

  %Q = [repdiag(L,dim) BT';BT CC];

  energy = 'spokes-and-rims';
  U = [];
  R = [];
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Energy','EvaluationPoint','Rotations'},{'energy','U','R'});
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

  % Rewrite energy in [Huang et al. 09] for arbitrary edge-sets
  % E(U) = ∑_k∈[1,r] ∑_{i,j}∈e_k w_ijk ‖ c_k × (p_i - p_j) - (u_i - u_j) ‖²

  % 3x3 blocks of crossproducts as matrices
  % Bik = [ ∑_j s_ijk w_ijk (p_i - p_j) ]×  where s_ijk is +1,-1 depending on
  % order of {ij} in e_k

  % 3x3 blocks of covariance matrices (sums of outerproducts)
  % **This is probably wrong in [Huang et al. 2009]**
  % Ckk = ∑_{ij}∈e_k w_ijk (p_i - p_j)(p_i - p_j)' 
  % 
  % Should probably also contain:
  % Ckk += ∑_{ij}∈e_k wijk ‖pi-pj‖² I

  % Number of vertices
  n = size(V,1);
  dim = size(V,2);
  % Simplex size
  ss = size(F,2);
  % Number of simplices
  m = size(F,1);

  L = -3*cotmatrix(V,F);

  % A(i,j) = w~=0 means that rotation edge set i contains edge allE(j,:) with
  % weight w
  LC = cotangent(V,F);

  % There seems to be an error (or at least omission) in [Huang et al. 2009].
  % The "C" matrix only computes part of each c × vij.
  % ‖Ck × Vij‖²
  % ‖Ck‖²‖Vij‖² - (Ck⋅Vij)²
  %  ^---diagonal part  ^------outer product (covariance matrix) part
  CSM = covariance_scatter_matrix(V,F,'Energy',energy);
  % Number of rotations
  nr = size(CSM,1)/dim;
  if isempty(U)
    % #R*dim by dim
    S = CSM*repmat(V,dim,1);
  else
    S = CSM*repmat(U,dim,1);
    SS = permute(reshape(S,[nr dim dim]),[2 3 1]);
    % S(:,:,r) := R(:,:,r) * S(:,:,r)
    if isempty(R)
      R = fit_rotations(SS);
    end
    RSS = zeros(size(R));
    for j = 1:dim
      RSS = RSS + bsxfun(@times,R(:,j,:),SS(j,:,:));
    end
    S = reshape(permute(RSS,[3 1 2]),[nr*dim dim]);
  end
  %S = reshape(permute(reshape(S,[],2,dim),[3 2 1]),dim,[])';

  X = cell(3,1);
  for d = 1:dim
    X{d} = arap_linear_block(V,F,d,'Energy',energy)';
  end
  Z = sparse(nr,n);
  BT = -3*[    Z  -X{3}  X{2}; ...
           X{3}     Z -X{1}; ...
          -X{2}  X{1}     Z];

  CI = repmat(1:nr*dim,dim,1)';
  CJ = repmat(reshape(1:nr*dim,nr,dim),dim,1);
  CC = -3*sparse(CI,CJ,S,nr*dim,nr*dim);
  CC = CC-diag(sparse(repmat(sum(reshape(diag(CC),nr,dim),2),dim,1)));

    % A X = B
    % But P' A P is easier to invert
    % P' A X = P' B
    % Let P Y = X --> Y = P' X
    % P' A P Y = P' B
    % Y = (P' A P) \ P' B
    % P' X  = (P' A P) \ P' B
    % X  = P (P' A P) \ P' B
    % inv(A) = P * inv(P' A P) P'

  %BCBT    = BT'*(CC\BT);
  P = sparse(1:nr*dim,reshape(1:nr*dim,dim,nr)',1);
  function X = chol_solve(A,B)
    [csL,p,csS] = chol(A);
    if p ~= 0
      warning('arap_hessian: chol failed');
      X = A \ B;
    else
      X = csS*(csL\(csL'\(csS'*B)));
    end
  end
  % TODO: is P'*CC*P always block diagonal?
  % yes, but extracting these blocks using a for loop is waaaay slower than
  % this:
  BCBT = BT'*(P*chol_solve(P'*CC*P,P'*BT));

  H = repdiag(L,dim) - BCBT;
end
