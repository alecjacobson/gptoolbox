function [E,ER] = arap_energy(V,F,U,varargin)
  % ARAP_ENERGY This function is meant as a very human readable groundtruth for
  % evaluating ARAP energies.
  % For faster energy computation see `arap_gradient.m`
  %
  % Inputs:
  %    V  #V by dim list of mesh vertex rest positions
  %    F  #F by simplex list of element indices into V
  %    U  #V by dim list of deformed vertex positions
  %   Optional:
  %     'Energy' followed by 'spokes','spokes-and-rims', or 'elements'
  %       {'elements' or 'spokes-and-rims'} for tets or triangles respectively.
  % Outputs:
  %    E  scalar energy
  %    ER  #R list of "per-rotation" energy contributions
  %
  %
  % Known issues:
  % If E_G is the energy returned by arap_gradient, then
  % 'spokes-and-rims' energy E = 3 E_G 
  % 'spokes' energy E = 2 E_G 
  % 'elements' energy E = E_G
  %


  switch size(F,2)
  case 4
    energy = 'elements';
  case 3
    energy = 'spokes-and-rims';
  end
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( {'Energy'}, ...
    {'energy'});
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

  dim = size(V,2);
  C = cotangent(V,F);

  n = size(V,1);
  m = size(F,1);
  switch energy
  case 'spokes'
    nr = size(V,1);
    edge_set = cell(nr,1);
    switch size(F,2)
    case 4
      AC = sparse( ...
        F(:,[2 3 1 4 4 4 3 1 2 1 2 3]),F(:,[3 1 2 1 2 3 2 3 1 4 4 4]), ...
        [C C],n,n);
    case 3
      AC = sparse(F(:,[2 3 1 3 1 2]),F(:,[3 1 2 2 3 1]),[C C],n,n);
    end
    % loop over vertices
    for r = 1:nr
      [Jr,~,Cr] = find(AC(:,r));
      Ir = repmat(r,numel(Jr),1);
      edge_set{r}.E = [Ir Jr];
      edge_set{r}.C = Cr;
    end
  case {'elements','spokes-and-rims'}
    nr = size(F,1);
    switch size(F,2)
    case 3
      I = F(:,[2 3 1]);
      J = F(:,[3 1 2]);
    case 4
      I = F(:,[2 3 1 4 4 4]);
      J = F(:,[3 1 2 1 2 3]);
    end
    edge_set = cell(nr,1);
    % loop over elements
    for r = 1:nr
      edge_set{r}.E = [I(r,:)' J(r,:)'];
      edge_set{r}.C = C(r,:);
    end
    if strcmp(energy,'spokes-and-rims')
      Fedge_set = edge_set;
      nr = size(V,1);
      edge_set = cell(nr,1);
      V2F = sparse( ...
        F,repmat(1:size(F,1),size(F,2),1)',repmat(1:size(F,2),size(F,1),1),n,m);
      % loop over vertices
      for r = 1:nr
        [~,Fr] = find(V2F(r,:));
        edge_set{r}.E = [];
        edge_set{r}.C = [];
        for f = Fr
          edge_set{r}.E = [edge_set{r}.E;Fedge_set{f}.E];
          edge_set{r}.C = [edge_set{r}.C Fedge_set{f}.C];
        end
      end

    end

  end

  ER = zeros(nr,1);
  % loop over rotations
  for r = 1:nr
    % get edges for this rotation
    Er = edge_set{r}.E;
    Cr = edge_set{r}.C;
    S = zeros(dim,dim);
    % loop over edges to build covariance matrix for rotation fitting
    for e = 1:size(Er,1)
      i = Er(e,1);
      j = Er(e,2);
      ev = V(j,:)-V(i,:);
      eu = U(j,:)-U(i,:);
      S = S + Cr(e) * eu'*ev;
    end
    R = fit_rotation(S);
    % Loop over edges to compute energy contribution
    for e = 1:size(Er,1)
      i = Er(e,1);
      j = Er(e,2);
      ev = V(j,:)-V(i,:);
      eu = U(j,:)-U(i,:);
      %ER(r) = ER(r) + 0.5 * Cr(e) * (eu - ev*R)*(eu - ev*R)';
      %ER(r) = ER(r) + 0.5 * Cr(e) * (eu*eu' - 2*ev*R*eu' + ev*R*R'*ev');
      ER(r) = ER(r) + 0.5 * Cr(e) * (eu*eu' - 2*ev*R*eu' + ev*ev');
    end
  end
  E = sum(ER);
end
