function [b,S] = bin_packing(a,V,varargin)
  % BIN_PACKING Given a list of n items of varying sizes (a), pack them into the
  % smallest number (b) of bins of equal size (V). 
  %
  % Inputs:
  %   a  #a list of item sizes
  %   V  scalar size of bins
  % Outputs:
  %   b  minimal number of bins
  %   S  #a list of indices into 1:b, assigning each item to a bin
  %
  % Example:
  %   a = ceil(rand(100,1)*10);
  %   V = 20;
  %   [b,S] = bin_packing(a,V);
  %   barh(full(sparse(S,1:numel(a),a)),'stacked');
  %   

  max_iters = 5;
  params_to_variables = containers.Map( ...
    {'MaxIters'},{'max_iters'});
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
  
  assert(all(a<V),'No item can be larger than single bin');
  n = numel(a);
  b = inf;
  S = zeros(n,1);
  % Try five times (first time is sort)
  for iter = 1:max_iters
    switch iter
    case 1
      P = (1:numel(a))';
    case 2
      [~,P] = sort(a,'descend');
    otherwise
      P = randperm(n);
    end
    ar = a(P);
    % start with no assignments
    Sr = zeros(n,1);
    % running remaining bin volumes
    Br = repmat(V,n,1);
    
    % naive O(n²) first fit:
    for i = 1:n
      % O(n) search to find first available bin
      j = find(Br>ar(i),1,'first');
      if j>=b
        break;
      end
      Sr(i) = j;
      Br(j) = Br(j) - ar(i);
    end
    if j >= b
      continue;
    end
    
    br = find(Br==V,1,'first')-1;
    if br < b
      b = br;
      B = Br;
      S(P) = Sr;
    end
  end


end
