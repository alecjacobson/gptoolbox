function [VV,FF,SS] = sqrt3(V,F,varargin)
  % SQRT3 Perform sqrt-3 subdivision. After n iterations of loop subivision, the
  % resulting mesh will have
  %   3^n |F| faces
  %   (3^(n+1)-3)/2 |F| + |E| edges
  %   (3^(n+1)-3)/6 |F| + |V| vertices
  %
  % [VV,FF,SS] = sqrt3(V,F)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by simplex-size list of simplex indices
  %   Optional:
  %     'Iterations' followed by number of recursive calls {1}
  % Outputs:
  %   VV  #VV by 3 list of output vertex positions
  %   FF  #FF by 3 list of triangle indices into VV
  %   SS  #VV by #V sparse matrix so that VV = SS*V

  % default values
  iters = 1;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Iterations'}, ...
    {'iters'});
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

  SS = speye(size(V,1),size(V,1));
  VV = V;
  FF = F;
  for iter = 1:iters
    n = size(SS,1);
    m = size(FF,1);
    C = n+(1:m)';
    E = edges(FF);
    A = adjacency_matrix(FF);
    [AI,AJ] = find(A);
    val = sum(A,1);
    An = (4-2*cos(2*pi./val))./9;
    Andval = An./val;
    A = sparse(AI,AJ,Andval(AI),n,n) + (diag(sparse(1-An)));
    S = [A;sparse(repmat(1:m,3,1)',FF,1/3,m,n)];
    SS = S*SS;
    FF = [FF(:,2) FF(:,3) C;FF(:,3) FF(:,1) C;FF(:,1) FF(:,2) C];
    %BC = barycenter(V,F);
    %VV = [V;BC];
    FF = flip_edges(FF,E,'NonConflicting',true);
  end
  VV = SS*V;

  % F₁ ← 3F₀
  % E₁ ← E₀ + 3F₀
  % V₁ ← V₀ + F₀

end
