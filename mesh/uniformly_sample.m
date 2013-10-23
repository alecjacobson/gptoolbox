function NV = uniformly_sample(V,un,varargin)
  % UNIFORMLY_SAMPLE  Uniformly sample an input curve 
  %
  % NV = uniformly_sample(V,un,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %  V  #V by dim list of ordered input points on an open curve
  %  un  # of samples {100}
  %  Optional:
  %    'Loop' followed by true or {false}
  %      NP = uniformly_sample(P,100,'Loop',true);
  %      % Equivalent to:
  %      NP = uniformly_sample(P([1:end 1],:),100);
  %      NP = NP(1:end-1,:);
  % Outputs:
  %  NV  un by dim list of ordered outputs points on an open curve
  %
  %

  loop = false;
  v = 1;
  while v<=numel(varargin)
    switch varargin{v}
    case 'Loop'
      v = v+1;
      assert(v<=numel(varargin));
      loop = varargin{v};
    otherwise
      error('Unsupported parameter "%s"',varargin{v});
    end
    v = v+1;
  end

  if loop
    V = V([1:end 1],:);
  end


  % Facets
  l = sqrt(sum((V(1:end-1,:)-V(2:end,:)).^2,2));
  totl = sum(l);

  % number of samples
  if nargin <2
    un = 100;
  end

  NV = zeros(un,2);
  v = 1;
  lv = 0;
  % TODO: There should be a vectorized way of doing this with cumsum
  for ni = 1:un-1
    li = ((totl-0)/(un-1))*(ni-1);
    while true
      if v<=size(l,1) && (lv+l(v))<=li
        lv = lv+l(v);
        v = v+1;
      else
        break;
      end
    end
    NV(ni,:) = V(v,:) + (li-lv)/l(v)*(V(mod(v,size(V,1))+1,:)-V(v,:));
  end
  NV(un,:) = V(end,:);

  if loop
    NV = NV(1:end-1,:);
  end
end
