function [SF,SVI,SV,A] = split_nonmanifold(F,varargin)
  % SPLIT_NONMANIFOLD Split a non-manifold (or non-orientable) mesh into a
  % manifold orientable mesh possibly with more connected components.
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


  % As if every corner became a vertex
  SF = reshape(1:numel(F),size(F));
  % These are edges between separated corners
  E = [SF(:,[2 3]);SF(:,[3 1]);SF(:,[1 2])];
  % F(E) is just the same as [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])]
  FE = F(E);
  FE
  % I(i) true indicates that the i-th edge in the original connected mesh has
  % found its consistently oriented twin
  [I,J] = ismember(FE,fliplr(FE),'rows');
  % E(I,:) Are edges between corners that have found their consistently oriented
  % twin in the original mesh

  % I think the idea is to connect faces which share compatibly oriented edges.
  % But this means for a non-manifold edge with 3 edges e,d,f we'd have e→d, e→f
  % but not d→f (w.l.o.g).

  % But I think ismember is just finding the 'first' occurrence.

  % This is just getting lucky!
  warning('this is bugged.');

  A = sparse(E(I,:),fliplr(E(J(I),:)),1,numel(F),numel(F));
  %[~,K] = conncomp(A);
  [~,K] = conncomp(A&A');

  [AI,AJ] = find(A);
  FE
  [AI AJ]
  fprintf('%d,%d → %d,%d\n',[FE(AI,:) FE(AJ,:)]');
  full(A)
  

  SVI = F(:);
  SV = V(SVI,:);

  SF = K(SF);
  SVI(K) = SVI;
  SV(K,:) = SV;

  [SV,IM,J] = remove_unreferenced(SV,SF);
  SF = J(SF);
  SVI = SVI(J);

  %if max(connected_components(SF)) == 1
  %keyboard
  %end

end
