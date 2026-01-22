function [NV,NI] = uniformly_sample(V,un,varargin)
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
  %    'E' followed by #E by 2 list of edges into V. In this case
  %    uniformly_sample will be run on each manifold patch. distributing un
  %    according to length of patch. Loop will be set accordingly
  % Outputs:
  %  NV  un by dim list of ordered outputs points on an open curve
  %  NI  un list of indices into NP of previous node w.r.t. arc length
  %
  %

  loop = false;
  v = 1;
  E = [];
  while v<=numel(varargin)
    switch varargin{v}
    case 'Loop'
      v = v+1;
      assert(v<=numel(varargin));
      loop = varargin{v};
    case 'E'
      v = v+1;
      assert(v<=numel(varargin));
      E = varargin{v};
    otherwise
      error('Unsupported parameter "%s"',varargin{v});
    end
    v = v+1;
  end

  if loop
    V = V([1:end 1],:);
  end

  if ~isempty(E)
    l = edge_lengths(V,E);
    totl = sum(l);
    [C] = manifold_patches(E);
    lC = accumarray(C(:),l);
    un_C = round(un*(lC./totl));
    % sum(un_C) not necessarily equal to un due to rounding
    S = sparse(1:size(E,1),C,1,size(E,1),max(C));
    NV = [];
    for c = 1:max(C)
      Ic = find(S(:,c));
      Ec = E(Ic,:);
      % Extract submesh. This next line is now O(|Ic|)
      [Vc,~,~,Ec] = remove_unreferenced(V,Ec,true);
      [IIc,J,K] = edges_to_path(Ec);
      loop_c = IIc(1) == IIc(end);
      [NVc,NIc] = uniformly_sample(Vc(IIc,:),un_C(c),'Loop',loop_c);
      % Not figuring these out right now.
      NIc(:) = nan;
      NV = [NV NVc'];
    end
    NV = NV';
    return;
  end


  % Facets
  l = sqrt(sum((V(1:end-1,:)-V(2:end,:)).^2,2));
  totl = sum(l);

  % number of samples
  if nargin <2
    un = 100;
  end

  NV = zeros(un,size(V,2));
  NI = zeros(size(NV,1),1);
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
    NI(ni) = v;
  end
  NV(un,:) = V(end,:);
  NI(un) = size(V,1);

  if loop
    NV = NV(1:end-1,:);
    NI = NI(1:end-1);
  end
end
