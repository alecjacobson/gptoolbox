function [G,E,R,data] = arap_gradient(V,F,U,varargin)
  % ARAP_GRADIENT Compute the Gradient of an as-rigid-as-possible for a mesh in
  % rest configuration (V,F) with vertices now at U, according to "A Simple
  % Geometric Model for Elastic Deformations" [Chao et al. 2010]
  %
  % G = arap_gradient(V,F,U)
  % [G,E,R] = arap_gradient(V,F,U,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by dim rest vertex positions
  %   F  #F by simplex-size simplex indices into V
  %   U  #V by dim deformed vertex positions
  %     Optional:
  %       'Rotations' followed by R  dim by dim by #R list of best fit
  %         rotations
  %       'Energy' followed by either 'spokes','spokes-and-rims','elements'
  % Outputs:
  %   G  #V by dim list of gradient vectors per vertex
  %   E  arap energy
  %   R  dim by dim by #R list of best fit rotations
  %   data 
  %
  % See also: arap, arap_hessian
  %
  % Example:
  % % Given a mesh (V,F) and deformed positions U0, flow to energy minimum
  % % using Newton's method.
  % clf;
  % hold on;
  %   tsurf(F,V,fphong,'FaceColor','r','SpecularStrength',0,'AmbientStrength',0.5);
  %   t = tsurf(F,U0,fphong,'FaceColor','b','SpecularStrength',0,'AmbientStrength',0.5);
  % hold off;
  % axis equal;
  % view(2);
  % camlight;
  % U = U0;
  % delta_t = 1e-1;
  % while true
  %   [G,E] = arap_gradient(V,F,U);
  %   U = U - delta_t * G;
  %   U = bsxfun(@plus,U,mean(V)-mean(U)+[max(V(:,1))-min(V(:,1)) 0 0]);
  %   t.Vertices = U;
  %   title(sprintf('E = %g\n',E),'FontSize',20);
  %   drawnow;
  % end
  %

  % default values
  switch size(F,2)
  case 4
    energy = 'elements';
  case 3
    energy = 'spokes-and-rims';
  end
  single_precision = true;
  R = [];
  data = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Energy','Rotations','Data','SinglePrecision'}, ...
    {'energy','R','data','single_precision'});
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

  % TODO: implement "flat" arap like `arap.m`, for now use placeholders with
  % non-flat defaults.
  ref_V = V;
  ref_F = F;
  dim = size(ref_V,2);
  flat = false;

  if isempty(data)
    ss = size(F,2);
    data.L = cotmatrix(V,F);
    if isempty(R)
      data.CSM = covariance_scatter_matrix(ref_V,ref_F,'Energy',energy);
    end
    [~,data.K] = arap_rhs(ref_V,ref_F,[],'Energy',energy);
  end

  % compute covariance matrix elements
  S = zeros(size(data.CSM,1),dim);
  S(:,1:dim) = data.CSM*repmat(U,dim,1);
  % dim by dim by n list of covariance matrices
  SS = permute(reshape(S,[size(data.CSM,1)/dim dim dim]),[2 3 1]);
  % fit rotations to each deformed vertex
  R = fit_rotations(SS,'SinglePrecision',single_precision);

  nr = size(R,3);
  Rcol = reshape(permute(R,[3 1 2]),nr*dim*dim,1);
  dV = data.K * Rcol;
  dV = reshape(dV,[size(V,1) dim]);

  G = -(data.L*U + dV);
  if nargout > 1
    E = trace(-U'*0.5*data.L*U - U'*dV - V'*0.5*data.L*V);
  end

end
