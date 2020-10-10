function [K,C,strain,A,M] = linear_elasticity_stiffness(V,F,varargin)
  %
  % [K,C,strain,A,M] = linear_elasticity_stiffness(V,F)
  %
  % Inputs:
  %   V  #V by d list of vertex positions
  %   F  #F by d+1 list of  element indices into V
  %   Optional:
  %     'Lambda'  followed by first Lamé parameter {1.7423333}, scalar
  %       (homogeneous) or #F by 1 list of per-element values
  %     'Mu'  followed by shear modulus {0.0115}, scalar (homogeneous) or #F by
  %       1 list of per-element values
  %     'Young'  followed by Young's modulus, scalar (homogeneous) or #F by 1
  %       list of per-element values
  %     'Nu'  followed by Poisson's ratio, scalar (homogeneous) or #F by 1 list
  %       of per-element values
  % Outputs:
  %   K  #V*d by #V*d sparse stiffness matrix
  %   C  #F**(d*(d+1)/2) by #F**(d*(d+1)/2) sparse constituitive model matrix 
  %   strain  #F*(d*(d+1)/2) by #V*d sparse strain matrix
  %   A  #F*(d*(d+1)/2) by #F*(d*(d+1)/2) diagonal element area matrix
  %   M  #V*d by #V*d sparse mass matrix
  %

  % Silicone rubber: http://www.azom.com/properties.aspx?ArticleID=920
  mu = 0.0115;
  % Bulk modulus
  K = 1.75;
  lambda = K-2/3*mu;
  young = [];
  nu = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Lambda','Mu','Nu','Young'}, ...
    {'lambda','mu','nu','young'});
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

  assert( ...
    (~isempty(lambda) && ~isempty(mu))||(~isempty(young) && ~isempty(nu)), ...
    'Must define either lambda and mu or young and nu');
  if (~isempty(young) && ~isempty(nu))
    lambda = young.*nu./((1+nu).*(1-2.*nu));
    mu = .5.*young./(1+nu);
  end

  dim = size(V,2);

  % This matches the matlab code by Jonas Koko, in
  % "Vectorized Matlab Codes for Linear Two-Dimensional Elasticity"

  % Gradient/divergence operator
  G = grad(V,F);

  I = speye(size(F,1));
  Z = sparse(size(F,1),size(V,1));
  lambda = diag(sparse(lambda));
  mu = diag(sparse(mu));
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
  switch dim
  case 2
    G1 = G(1:size(F,1),:);
    G2 = G(size(F,1)+(1:size(F,1)),:);
    % 3#F by 2#V
    strain = [G1 Z;Z G2;G2 G1];
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
    %C = kroneye(C,size(F,1));
    C = [ ...
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
    A = diag(sparse(doublearea(V,F)/2));
  case 3
    G1 = G(1:size(F,1),:);
    G2 = G(size(F,1)+(1:size(F,1)),:);
    G3 = G(2*size(F,1)+(1:size(F,1)),:);
    strain = [G1 Z Z;Z G2 Z; Z Z G3; Z G3 G2; G3 Z G1; G2 G1 Z];
    C = [ ...
        (lambda+2*mu)*I        lambda*I        lambda*I     0*I     0*I     0*I ; ...
        lambda*I (lambda+2*mu)*I        lambda*I     0*I     0*I     0*I ; ...
        lambda*I        lambda*I (lambda+2*mu)*I     0*I     0*I     0*I ; ...
        0*I             0*I             0*I    mu*I     0*I     0*I ; ...
        0*I             0*I             0*I     0*I    mu*I     0*I ; ...
        0*I             0*I             0*I     0*I     0*I    mu*I ];
    A = diag(sparse(volume(V,F)));
  end
  Z = sparse(size(V,1),size(F,1));
  D = strain';
  A = repdiag(A,dim*(dim+1)/2);
  K = D * A * C * strain;


  M = massmatrix(V,F);
  M = repdiag(M,size(V,2));

end
