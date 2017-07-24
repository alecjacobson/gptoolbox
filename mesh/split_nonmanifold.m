function [SF,SVI,SV] = split_nonmanifold(F,varargin)
  % SPLIT_NONMANIFOLD Split a non-manifold mesh into a manifold mesh possibly
  % with more connected components.
  %
  % Inputs:
  %   F  #F by 3 list of input tringle indices into some vertex list V
  %   Optional:
  %     'V' followed by #V by dim vertex list
  % Outputs:
  %   SF  #F by 3 list of output tringle indices into V(SVI,:)
  %   SVI  #SV list of indices into V identifying vertex positions
  %   SV  #SV by dim list of vertex positions so that SV = V(SVI,:)
  %

  V = [];
  params_to_variables = containers.Map( ...
    {'V'}, ...
    {'V'});
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

  if isempty(V)
    V = zeros(max(F(:)),0);
  end


  SF = reshape(1:numel(F),size(F));
  E = [SF(:,[2 3]);SF(:,[3 1]);SF(:,[1 2])];
  [I,J] = ismember(F(E),fliplr(F(E)),'rows');
  A = sparse(E(I,:),fliplr(E(J(I),:)),1,numel(F),numel(F));
  [~,K] = conncomp(A);

  SVI = F(:);
  % SV = V(VI,:);
  SV = V(SVI,:);

  BC = barycenter(SV,SF);
  SF = K(SF);
  SVI(K) = SVI;
  SV(K,:) = SV;

  [SV,IM,J] = remove_unreferenced(SV,SF);
  SF = J(SF);
  SVI = SVI(J);

end
