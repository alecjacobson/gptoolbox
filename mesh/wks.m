function [K,data] = wks(V,F,varargin)
  % WKS Wave Kernel Signature "TheWave Kernel Signature: A Quantum Mechanical
  % Approach to Shape Analysis" [Aubry et al. 2011]
  %
  % Inputs:
  %   V  #V by dim list of mesh vertex positions
  %   F  #F by ss list of element indices into V
  % Outputs:
  %   K  #V by #K list of HKS feature vectors for each vertex
  %   data
  %     E  #E list of eigen values
  %     phi  #V by #E list of eigen functions
  %     phi_sqr  #V by #E list of squared eigen functions
  %     Ce  #E list of exp(- eigen values)
  %     M  #V by #V mass matrix
  %     L  #V by #V Laplacian matrix
  %     e  #e list of energy values
  % 
  % Example:
  %   K = wks(V,F);
  %   % Compare all vertices' signatures to first using metric of [Aubry et al.
  %   % 2011]
  %   tsurf(F,V,'CData',sum(abs((K-K(1,:))./(K+K(1,:))),2),fphong);;

  emin = [];
  emax = [];
  k = 300;
  nt = 100;
  sigma = [];
  data = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'emax','emin','Data','sigma','NumModes','NumTimeSteps'}, ...
    {'emax','emin','data','sigma','k','nt'});
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
    [data.phi,data.E] = eigs(-data.L,data.M,k,'sm');
    data.E = diag(data.E);
    data.phi_sqr = data.phi.^2;
  end
  % The parameter section of WKS paper is confusing: "In all our experiments,
  % the parameters were fixed..."
  %
  % This is implies that sigma, emin, and emax are the solution to a linear
  % system:
  log_E_min = log(min(data.E(data.E>1e-8)));
  log_E_max = log(max(data.E(data.E>1e-8)));
  %if isempty(sigma)
  %  delta =  (log_E_max-log_E_min)/(nt-1+28);
  %  sigma = 7*delta;
  %end
  %if isempty(emin)
  %  emin = log_E_min-2*sigma;
  %end
  %if isempty(emax)
  %  emax = log_E_max+2*sigma;
  %end

  assert(isempty(emin));
  assert(isempty(emax));
  assert(isempty(sigma));
  sol = [ ...
    1/(nt-1) -1/(nt-1)  0 -1; ...
    0        0         -1  7; ...
    1        0          2  0; ...
    0        1         -2  0] \ ...
    [0;0;log_E_max;log_E_min];
  emax = sol(1);
  emin = sol(2);
  sigma = sol(3);
  delta = sol(4);

  data.e = linspace(emin,emax,nt);

  data.fE = exp((-(data.e-log(data.E)).^2)/(2*sigma^2));
  data.Ce = sum(data.fE,1).^-1;
  K = (data.phi_sqr * data.fE) ./ data.Ce;

end
