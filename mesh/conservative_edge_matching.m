function [EE,I] = conservative_edge_matching(E_orig,varargin)
  % CONSERVATIVE_EDGE_MATCHING  Find set of given edges so that no two edges
  % share a vertex.
  %
  % [EE,I] = conservative_edge_matching(E)
  % 
  % Inputs:
  %   E  #E by 2 list of edges 
  %   Optional:
  %      'Method'  followed by one of the following:
  %        'one-per-component'  simply take one edge per conneceted component
  %          of edges
  %        'random'  use a randomized algorithm to try to grab many per
  %          component
  % Outputs:
  %   EE  #EE by 2 list of edges
  %   I  #EE list of indices into E, so that EE = E(I,:)
  %

  if size(E_orig,1) == 0
    EE = E_orig;
    I = [];
    C = [];
    return;
  end


  % default values
  method = 'one-per-component';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Method'}, ...
    {'method'});
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

  % Do this in reduced graph of just vertices incident on E. Otherwise this is
  % O(max(E(:))) rather than O(size(E))
  [~,~,E] = unique(E_orig(:));
  E = reshape(E,size(E_orig));


  switch method
  case 'one-per-component'
    % Find one edge per connected component of edges (ideally we would find a
    % "perfect matching")
    ne = size(E,1);
    C = connected_components(E);
    C = C(E(:,1));
    % Use reverse order so we can take max
    E2C = sparse(1:ne,C,ne-(1:ne)+1,ne,max(C));
    I = max(E2C,[],1);
    I = I(I>0);
    % reverse
    I = ne - I + 1;
  case 'recursive'
    keep = 1:size(E,1);
    Ekeep = E(keep,:);
    J = 1:size(E,1);
    I = [];
    while true
      % This should use random for large meshes
      [~,Ikeep] = conservative_edge_matching(Ekeep,'Method','random');
      I = [I J(Ikeep)];
      % Remove edges connected to E
      keep = find(~any(ismember(Ekeep,Ekeep(Ikeep,:)),2));
      if isempty(keep)
        break;
      end
      % Indices in original edges
      J = J(keep);
      Ekeep = Ekeep(keep,:);
    end

  case 'random'
    ns = size(E,1);
    while true
      I = randperm(size(E,1));
      I = I(1:ns);
      % Reduced graph of randomly selected edges
      [~,~,EI] = unique(E(I,:));
      % Check if every vertex is valence 1
      V = full(sparse(EI,1,1));
      if max(V) == 1
        break;
      end
      ns = floor(ns/2);
      if ns <= 1
        I = 1;
        break;
      end
    end
  end

  % Output on original graph
  EE = E_orig(I,:);
end
