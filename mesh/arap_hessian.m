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
  % 
  %Q = [repdiag(L,dim) BT';BT CC];

  energy = 'spokes-and-rims';
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'Energy'},{'energy'});
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

  switch ss
  case 3
    L = -3*cotmatrix(V,F);
  case 4
    L = -3*cotmatrix3(V,F);
  end

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
  S = CSM*repmat(V,dim,1);
  SS = permute(reshape(S,[size(CSM,1)/dim dim dim]),[2 3 1]);
  %S = reshape(permute(reshape(S,[],2,dim),[3 2 1]),dim,[])';
  CI = repmat(1:nr*dim,dim,1)';
  CJ = repmat(reshape(1:nr*dim,nr,dim),dim,1);
  CC = -3*sparse(CI,CJ,S,nr*dim,nr*dim);
  CC = CC-diag(sparse(repmat(sum(reshape(diag(CC),nr,dim),2),dim,1)));

  X = cell(3,1);
  for d = 1:dim
    X{d} = arap_linear_block(V,F,d,'Energy',energy)';
  end
  Z = sparse(nr,n);
  BT = -3*[    Z  -X{3}  X{2}; ...
           X{3}     Z -X{1}; ...
          -X{2}  X{1}     Z];

  % Inverting CC explicitly is also possible (and fast if done correctly) but
  % `chol` seems to be pretty fast.
  [LCC,p]=chol(CC);
  BCBT = BT'*(LCC\(LCC'\BT));
  H = repdiag(L,dim) - BCBT;
end
