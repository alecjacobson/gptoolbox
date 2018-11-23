function [K,MK,data] = hks(V,F,varargin)
  % HKS Heat Kernel Signature "A Concise and Provably Informative Multi-Scale
  % Signature Based on Heat Diffusion" [Sun et al. 2009]
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by ss list of element indices into V
  % Outputs:
  %   K  #V by #K list of HKS feature vectors for each vertex
  %   MK  #K list of mass-weighted spatial sums
  %   data
  %     lambda  #lambda list of eigen values
  %     phi  #V by #lambda list of eigen functions
  %     phi_sqr  #V by #lambda list of squared eigen functions
  %     exp_m_lambda  #lambda list of exp(- eigen values)
  %     M  #V by #V mass matrix
  %     L  #V by #V Laplacian matrix
  % 
  % Example:
  %   [K,MK] = hks(V,F);
  %   % Compare all vertices' signatures to first using metric of [Sun et al.
  %   % 2009]
  %   tsurf(F,V,'CData',sqrt(sum(((K-K(1,:))./MK).^2,2)),fphong);

  tmin = [];
  tmax = [];
  data.t = [];
  data.exp_m_lambda = [];
  data.phi_sqr = [];
  data.M = [];
  k = 300;
  nt = 100;
  data = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'tmax','tmin','Data','NumModes','NumTimeSteps'}, ...
    {'tmax','tmin','data','k','nt'});
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

  if isempty(data)
    data.L = cotmatrix(V,F);
    data.M = massmatrix(V,F);
    [data.phi,data.lambda] = eigs(-data.L,data.M,k,'sm');
    data.lambda = diag(data.lambda);
    data.phi_sqr = data.phi.^2;
    data.exp_m_lambda = exp(-data.lambda);
  end
  if isempty(tmin)
    tmin = 4*log(10)/max(data.lambda);
  end
  if isempty(tmax)
    tmax = 4*log(10)/min(data.lambda(data.lambda>min(data.lambda)));
  end

  data.t = exp(linspace(log(tmin),log(tmax),nt));

  K = (data.phi_sqr)*(data.exp_m_lambda.^data.t);
  MK = diag(data.M)'*K;

end
