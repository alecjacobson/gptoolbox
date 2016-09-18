function [U,Ud,data] = linear_elasticity(V,F,b,bc,varargin)
  % LINEAR_ELASTICITY Compute the deformation of a 2D solid object according to
  % a linear model of elasticity, assuming a linear isotropic material.
  % 
  % U = linear_elasticity(V,F,b,bc)
  % [U,Ud,K] = linear_elasticity(V,F,b,bc,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 2 list of vertex positions
  %   F  #F by 3 list of triangle element indices into V
  %   b  #b list of indices into V of fixed vertices
  %   bc #bc by 2 list of fixed vertex positions
  %   Optional:
  %     'Lambda'  followed by first Lamé parameter {1.7423333}, scalar
  %       (homogeneous) or #F by 1 list of per-element values
  %     'Mu'  followed by shear modulus {0.0115}, scalar (homogeneous) or #F by
  %       1 list of per-element values
  %     'Young'  followed by Young's modulus, scalar (homogeneous) or #F by 1
  %       list of per-element values
  %     'Nu'  followed by Poisson's ratio, scalar (homogeneous) or #F by 1 list
  %       of per-element values
  %     'U0'  followed by #V by 2 list of previous displacements
  %     'Ud0'  followed by #V by 2 list of previous velocities: (U0 - Um1)/dt
  %     'BodyForces'  followed by #V by 2 list of body forces
  %     'TimeStep' followed by time step {0.1}
  % Outputs:
  %   U  #V by 2 list of vertex displacements
  %   Ud  #V by 2 list of vertex velocities
  %   data  precomputation data
  %     data.A  #F*3 by #F*3 diagonal element area matrix
  %     data.K  #V*2 by #V*2 sparse stiffness matrix
  %     data.M  #V*2 by #V*2 sparse mass matrix
  %     data.strain  #F*3 by #V*2 sparse strain matrix
  %     data.dt  timestep
  %     data.C  #F*3 by #F*3 sparse constituitive model matrix 
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

  assert(size(V,2) == 2,'Only 2D meshes are supported');

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

  if isempty(data)
    tic;
    data.dt = dt;
    assert( ...
      (~isempty(lambda) && ~isempty(mu))||(~isempty(young) && ~isempty(nu)), ...
      'Must define either lambda and mu or young and nu');
    if (~isempty(young) && ~isempty(nu))
      lambda = young.*nu./((1+nu).*(1-2.*nu));
      mu = .5.*young./(1+nu);
    end

    % This matches the matlab code by Jonas Koko, in
    % "Vectorized Matlab Codes for Linear Two-Dimensional Elasticity"

    % Gradient/divergence operator
    G = grad(V,F);

    % Strain tensor 
    %
    %   ϵ = ½(∇u + (∇u)')
    %   ϵ = ½ // ∂u₁/∂x₁  ∂u₂/∂x₁ \  + / ∂u₁/∂x₁  ∂u₁/∂x₂ \\
    %         \\ ∂u₁/∂x₂  ∂u₂/∂x₂ /    \ ∂u₂/∂x₁  ∂u₂/∂x₂ //
    %
    %                                "Voigt" notation
    %   ϵ₁₁ = ∂u₁/∂x₁              = ϵ₁
    %   ϵ₂₂ = ∂u₂/∂x₂              = ϵ₂
    %   ϵ₁₂ = ½(∂u₂/∂x₁ + ∂u₁/∂x₁) = ½ ϵ₃
    %   ϵ₂₁ = ϵ₁₂                  = ½ ϵ₃
    %  
    G1 = G(1:size(F,1),:);
    G2 = G(size(F,1)+(1:size(F,1)),:);
    Z = sparse(size(F,1),size(V,1));
    % 3#F by 2#V
    data.strain = [G1 Z;Z G2;G2 G1];

    % Stiffness tensor
    %
    %    σ = C:ϵ        %  A:B = Aij Bij 
    %                   %      = ∑∑ Aij Bij, where in this case Aij is a 2x2
    %                   %                    matrix, and Bij is a scalar
    %  
    % For each face we have:
    %   
    %    2x2 = 2x2x2x2 2x2
    %    σf  = Cf : ϵf
    %    σ = ∑∑ Cij ϵij, where Cij is a 2x2 matrix
    %    σkl = ∑∑ Cijkl ϵij, where Cijkl is a scalar
    %
    % But really ϵf and σf are just 3 distinct values:
    %
    %   σ₁ = [ϵ₁ ϵ₂ ϵ₃] [ c₁₁ ; c₁₂ ; c₁₃ ]
    %   σ₂ = [ϵ₁ ϵ₂ ϵ₃] [ c₂₁ ; c₂₂ ; c₂₃ ]
    %   σ₃ = [ϵ₁ ϵ₂ ϵ₃] [ c₃₁ ; c₃₂ ; c₃₃ ]
    %
    %    /σ₁\     /c₁₁ c₁₂ c₁₃\  /ϵ₁\
    %   | σ₂ | = | c₂₁ c₂₂ c₂₃ || ϵ₂ |
    %    \σ₃/     \c₃₁ c₃₂ c₃₃/  \ϵ₃/
    %  
    % So if σ is a 3#F by 1 vector and ϵ is a 3#F vector then:
    %  
    %   σ = C ϵ
    %        /C₁₁ C₁₂ C₁₃\  /ϵ₁\
    %   σ = | C₂₁ C₂₂ C₂₃ || ϵ₂ |
    %        \C₃₁ C₃₂ C₃₃/  \ϵ₃/
    % 
    %  where C is 3#F by 3#F matrix and Cij = diagonal #F by #F matrix.
    %
    % For Isotropic homogeneous media, we have that:
    %
    %   σij = λ δij ϵkk + 2μ ϵij
    %   σij = λ δij (∑ ϵkk) + 2μ ϵij
    % 
    % where λ is Lamé's first parameter and μ is the shear modulus: the bulk
    % modulus is thus K := λ + ⅔ μ
    %
    % Or in Voigt notation:
    % 
    %   σ₁ = σ₁₁ = λ (ϵ₁ + ϵ₂) + 2μ ϵ₁
    %   σ₂ = σ₂₂ = λ (ϵ₁ + ϵ₂) + 2μ ϵ₂
    %   σ₃ = σ₁₂ = λ (ϵ₁ + ϵ₂) + 2μ ϵ₁₂
    %            = λ (ϵ₁ + ϵ₂) + μ ϵ₃
    %
    %        //λ  λ 0\   /2μ  0  0\\  
    %  σ =  || λ  λ 0 |+|  0 2μ  0 || ϵ
    %        \\0  0 0/   \ 0  0  μ//
    %
    %
    %Z = sparse(size(F,1),size(F,1));
    %I = speye(size(F,1));
    %C = lambda*[[I I Z;I I Z;Z Z Z]] + mu*[2*I Z Z;Z 2*I Z;Z Z I];
    %C = lambda*[1 1 0;1 1 0;0 0 0] + mu*diag([2 2 1]);
    %data.C = kroneye(C,size(F,1));
    I = speye(size(F,1));
    lambda = diag(sparse(lambda));
    mu = diag(sparse(mu));
    data.C = [ ...
      (lambda+2*mu)*I        lambda*I  0*I; ...
             lambda*I (lambda+2*mu)*I  0*I; ...
                  0*I             0*I mu*I];

    %   ∇⋅σ = /∇⋅/σ₁₁\  ∇⋅/σ₁₂\\
    %         \  \σ₂₁/    \σ₂₂//
    % 
    % If D is the divergence operator then D is 2#V by 3#F, where σ is 3#F by 1
    % vectorized stress tensor using Voigt notation:
    %
    %   X = D σ
    %
    Z = sparse(size(V,1),size(F,1));
    D = [G1' Z G2';Z G2' G1'];
    A = diag(sparse(doublearea(V,F)/2));
    data.A = blkdiag(A,A,A);
    data.K = D * data.A * data.C * data.strain;

    data.M = massmatrix(V,F);
    data.M = repdiag(data.M,size(V,2));

    % ∇⋅σ + F = ü
    % Ku₂ + MF = M(u₂-2u₁+u₀)/dt²
    % dt²Ku₂ + dt²MF = M(u₂-2u₁+u₀)
    % (dt²K - M)u₂  = -dt²MF + M(-2u₁+u₀)
    % -(dt²K - M)u₂  = dt²MF - M(-2u₁+u₀)
    % (M-dt²K)u₂  = dt²MF + M(2u₁-u₀)
    % (M-dt²K)u₂  = M*(dt²F + 2u₁-u₀)
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


    % B = data.M*(data.dt^2*fext(:) + 2*U0(:) - Um1(:));
    B = data.M*(data.dt^2*fext(:) + U0(:) + data.dt*Ud0(:));
    % Fix each coordinate
    bb = reshape(bsxfun(@plus,reshape(b,[],1),(0:size(V,2)-1)*size(V,1)),1,[]);
    % solve once to set data.mqwf
    [U,data.mqwf] = min_quad_with_fixed(A,-2*B,bb,bc(:));
    data.solve  = @(B,bc) min_quad_with_fixed(A,-2*B,bb,bc(:),[],[],data.mqwf);
  else
    B = data.M*(data.dt^2*fext(:) + U0(:) + data.dt*Ud0(:));
    U = data.solve(B,bc(:));
  end
  U = reshape(U,size(V));
  Ud = (U-U0)/data.dt;

end
