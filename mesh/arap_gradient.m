function G = arap_gradient(V,F,U,varargin)
  % ARAP_GRADIENT Compute the Gradient of an as-rigid-as-possible for a mesh in
  % rest configuration (V,F) with vertices now at U, according to "A Simple
  % Geometric Model for Elastic Deformations" [Chao et al. 2010]
  %
  % Inputs:
  %   V  #V by dim rest vertex positions
  %   F  #F by simplex-size simplex indices into V
  %   U  #V by dim deformed vertex positions
  % Outputs:
  %   G  #V by dim list of gradient vectors per vertex
  %

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

  ss = size(F,2);
  switch ss
  case 3
    L = cotmatrix(V,F);
  case 4
    L = cotmatrix3(V,F);
  otherwise
    error('F must be #F by 3 or #F by 4');
  end

  % TODO: implement "flat" arap like `arap.m`, for now use placeholders with
  % non-flat defaults.
  ref_V = V;
  ref_F = F;
  dim = size(ref_V,2);
  flat = false;

  CSM = covariance_scatter_matrix(ref_V,ref_F,'Energy',energy);
  % compute covariance matrix elements
  S = zeros(size(CSM,1),dim);
  S(:,1:dim) = CSM*repmat(U,dim,1);
  % dim by dim by n list of covariance matrices
  SS = permute(reshape(S,[size(CSM,1)/dim dim dim]),[2 3 1]);
  % fit rotations to each deformed vertex
  R = fit_rotations(SS);
  dV = arap_rhs(ref_V,ref_F,R,'Energy',energy);

  G = -(L*U + dV);

end
