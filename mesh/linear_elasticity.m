function [U,Ud,data] = linear_elasticity(V,F,b,bc,varargin)
  % LINEAR_ELASTICITY Compute the deformation of a 2D solid object according to
  % a linear model of elasticity, assuming a linear isotropic material.
  % 
  % U = linear_elasticity(V,F,b,bc)
  % [U,Ud,K] = linear_elasticity(V,F,b,bc,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by d list of vertex positions
  %   F  #F by d+1 list of element indices into V
  %   b  #b list of indices into V of fixed vertices
  %   bc #bc by d list of fixed vertex positions
  %   Optional:
  %     'Lambda'  followed by first Lamé parameter {1.7423333}, scalar
  %       (homogeneous) or #F by 1 list of per-element values
  %     'Mu'  followed by shear modulus {0.0115}, scalar (homogeneous) or #F by
  %       1 list of per-element values
  %     'Young'  followed by Young's modulus, scalar (homogeneous) or #F by 1
  %       list of per-element values
  %     'Nu'  followed by Poisson's ratio, scalar (homogeneous) or #F by 1 list
  %       of per-element values
  %     'U0'  followed by #V by d list of previous displacements
  %     'Ud0'  followed by #V by d list of previous velocities: (U0 - Um1)/dt
  %     'BodyForces'  followed by #V by d list of body forces
  %     'TimeStep' followed by time step {0.1}
  %     'Data'  see output {[]}
  % Outputs:
  %   U  #V by d list of vertex displacements
  %   Ud  #V by d list of vertex velocities
  %   data  precomputation data
  %     data.A  #F*(d*(d+1)/2) by #F*(d*(d+1)/2) diagonal element area matrix
  %     data.K  #V*d by #V*d sparse stiffness matrix
  %     data.M  #V*d by #V*d sparse mass matrix
  %     data.strain  #F*(d*(d+1)/2) by #V*d sparse strain matrix
  %     data.dt  timestep
  %     data.C  #F**(d*(d+1)/2) by #F**(d*(d+1)/2) sparse constituitive model matrix 
  %     data.mqwf  precomputation for implicit solve (from min_quad_with_fixed)
  %     data.solve  function handle for conducting implicit step
  % 
  % Example:
  %   % Fit to half the unit square
  %   V = V/(2*max(max(V)-min(V))); 
  %   % Initialize as stretched object
  %   U = 1.5*V-V;
  %   Ud = zeros(size(V));
  %   t = tsurf(F,V+U);
  %   axis equal;
  %   axis manual;
  %   while true
  %     [U,Ud] = linear_elasticity(V,F,[],[],'U0',U,'Ud0',Ud);
  %     t.Vertices = V+U;
  %     drawnow;
  %   end
  %   

  data = [];
  % Time step
  dt = 0.1;
  % Silicone rubber: http://www.azom.com/properties.aspx?ArticleID=920
  mu = 0.0115;
  % Bulk modulus
  K = 1.75;
  lambda = K-2/3*mu;
  young = [];
  nu = [];
  U0 = zeros(size(V));
  Ud0 = zeros(size(V));
  fext = zeros(size(V));
  %% Parameters so that off-diagonals _should_ be zero
  %lambda = 1;
  %mu = -2*lambda;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Lambda','Mu','Nu','Young','U0','Ud0','BodyForces','TimeStep','Data'}, ...
    {'lambda','mu','nu','young','U0','Ud0',      'fext',      'dt','data'});
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


  first_solve = false;
  if isempty(data)
    first_solve = true;
    tic;
    data.dt = dt;
    assert( ...
      (~isempty(lambda) && ~isempty(mu))||(~isempty(young) && ~isempty(nu)), ...
      'Must define either lambda and mu or young and nu');
    if (~isempty(young) && ~isempty(nu))
      lambda = young.*nu./((1+nu).*(1-2.*nu));
      mu = .5.*young./(1+nu);
    end
    % young = mu*(3*lambda+2*mu)/(lambda+mu);
    % nu = lambda/(2*(lambda+mu));

    [data.K,data.C,data.strain,data.A,data.M] = ...
      linear_elasticity_stiffness(V,F,'Lambda',lambda,'Mu',mu);


    % ∇⋅σ + F = ü
    % Ku₂ + MF = M(u₂-2u₁+u₀)/dt²
    % dt²Ku₂ + dt²MF = M(u₂-2u₁+u₀)
    % (dt²K - M)u₂  = -dt²MF + M(-2u₁+u₀)
    % -(dt²K - M)u₂  = dt²MF - M(-2u₁+u₀)
    % (M-dt²K)u₂  = dt²MF + M(2u₁-u₀)
    % (M-dt²K)u₂  = M*(dt²F + 2u₁-u₀)
    % ud₀ = (u₁-u₀)/dt
    % dt*ud₀ = u₁-u₀
    % (M-dt²K)u₂  = M*(dt²F + u₁ + u₁-u₀)
    % (M-dt²K)u₂  = M*(dt²F + u₁ + dt*ud₀)
    A = data.M+data.dt^2*data.K;
    % ud = (u - u0)/dt
    % udd = ((u - u0)-(u0-um1)/dt²
    % udd = (u - 2u0 +um1)/dt²
    % ud*dt = (u - u0)
    % ud*dt - u = -u0
    % u - ud*dt = u0
    % udd*dt² = u - 2u0 +um1
    % udd*dt² - u + 2u0 = um1
    % 2u0-um1
    % 2(u - ud*dt)-(udd*dt² - u + 2u0)
    % 2u - 2ud*dt-udd*dt² + u - 2u0
    % 3u - 2ud*dt-udd*dt² - 2(u - ud*dt)
    % 3u - 2ud*dt-udd*dt² - 2u + ud*dt
    % u-ud*dt-udd*dt²


  end
  B = data.M*(data.dt^2*fext(:) + U0(:) + data.dt*Ud0(:));
  if first_solve
    % Fix each coordinate
    bb = reshape(bsxfun(@plus,reshape(b,[],1),(0:size(V,2)-1)*size(V,1)),1,[]);
    % solve once to set data.mqwf
    [U,data.mqwf] = min_quad_with_fixed(A,-2*B,bb,bc(:));
    data.solve  = @(B,bc) min_quad_with_fixed(A,-2*B,bb,bc(:),[],[],data.mqwf);
  else
    U = data.solve(B,bc(:));
  end
  U = reshape(U,size(V));
  Ud = (U-U0)/data.dt;

end
