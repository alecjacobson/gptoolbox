function [N,I,B,r] = random_points_on_mesh(V,F,n,varargin)
  % RANDOM_POINTS_ON_MESH Uniform random sampling of a mesh.
  %
  % N = random_points_on_mesh(V,F,n)
  % [N,I,B,r] = random_points_on_mesh(V,F,n,'ParameterName',parameter_value,...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   n  number of samples
  %   Optional:
  %     'Color' followed by "color of noise", one of:
  %       {'white'}  uniform random noise
  %       'blue'  approximate blue noise via approximate Poisson disk sampling
  %     'MaxIter'  followed by maximum number of iterations, only applicable
  %       for blue noise.
  %     'Bias'  **experimental** followed by #F list of biases per triangle
  %       {doublearea(V,F)}.  This bias is a little funny because it effects
  %       the sampling, but not the distances...probably it gets "washed out"
  %       by the iterations.
  % Outputs:
  %   N  n by dim list of sample positions
  %   I  n list of face indices upon which N are
  %   B  n by 3 list of barycentric coordinates so that:
  %     N =  ...
  %       bsxfun(@times,B(:,1),V(F(I,1),:)) +  ...
  %       bsxfun(@times,B(:,2),V(F(I,2),:)) +  ...
  %       bsxfun(@times,B(:,3),V(F(I,3),:));
  %   r  radius used for blue noise sampling
  %
  % See also: farthest_points
  %

  % default values
  color = 'white';
  max_iter = 100;
  A = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Color','MaxIter','Bias'}, ...
    {'color','max_iter','A'});
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

  if isempty(A)
    A = doublearea(V,F);
  end

  assert(size(A,1) == size(F,1));

  switch color
  case 'white'
    A = A./sum(A);
    A0 = [0;A];
    C = cumsum(A0);
    R = rand(n,1);
    [~,I] = histc(R,C);
    S = rand(n,1);
    T = rand(n,1);
    B = [1-sqrt(T) (1-S).*sqrt(T) S.*sqrt(T)];
    N =  ...
      bsxfun(@times,B(:,1),V(F(I,1),:)) +  ...
      bsxfun(@times,B(:,2),V(F(I,2),:)) +  ...
      bsxfun(@times,B(:,3),V(F(I,3),:));
  case 'blue'
    % This is a "cheap hack" way of getting something like Poisson-Disk
    % sampling which approximates a blue noise sampling.
    k = 10;
    %BV = biharmonic_embedding(V,F,k);
    randsample = @(n) random_points_on_mesh(V,F,n,'Color','white','Bias',A);
    embed = @(N,I,B) ...
      bsxfun(@times,B(:,1),V(F(I,1),:)) +  ...
      bsxfun(@times,B(:,2),V(F(I,2),:)) +  ...
      bsxfun(@times,B(:,3),V(F(I,3),:));      
    % initial sampling
    [N,I,B] = randsample(n);
    BN = embed(N,I,B);
    % Expected radius
    [~,D] = knnsearch(BN,BN,'K',2);
    D = D(:,2);
    %r = mean(D);
    r = mean(D)*sqrt(2);

    for iter = 1:max_iter
      % Maybe `rangesearch` would be faster
      %
      % Should only query new points...
      [kI,D] = knnsearch(BN,BN,'K',2);
      % Discard first match: always self
      kI = kI(:,2);
      D = D(:,2);
      % Find darts to rethrow, first of each pair with small distance
      rethrow = unique(kI( kI>(1:numel(kI))' & D<r));
      %fprintf('iter: %d, #rethrow: %d\n',iter,numel(rethrow));
      if numel(rethrow) == 0
        break;
      end
      [N(rethrow,:),I(rethrow),B(rethrow,:)] = randsample(numel(rethrow));
      keep = setdiff(1:n,rethrow);
      N = N([keep(:);rethrow(:)],:);
      I = I([keep(:);rethrow(:)],:);
      B = B([keep(:);rethrow(:)],:);
      BN = embed(N,I,B);
      if iter == max_iter
        warning('Max iteration (%d) reach without convergence',max_iter);
      end
    end

  end

end
