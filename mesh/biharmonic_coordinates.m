function [W,A,K,M,L,N] = biharmonic_coordinates(V,Ele,b,varargin)
  % BIHARMONIC_COORDINATES  Compute the linearly precise generalized
  % coordinates as described in "Linear Subspace Design for Real-Time Shape
  % Deformation" [Wang et al. 2015] **not** to be confused with "Bounded
  % Biharmonic Weights ..." [Jacobson et al. 2011] or "Biharmonic Coordinates"
  % [Weber et al. 12].
  %
  % W = biharmonic_coordinates(V,Ele,b)
  % [W,A,K,M,L,N] = biharmonic_coordinates(V,Ele,b)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   Ele  #Ele by dim+1 list of simplicial element indices into V
  %   b  #b list of indices into V of control points
  % Outputs:
  %   W  #V by #b list of coordinates
  %   A  #V by #V quadratic coefficients matrix
  %   K  #V by #V "Laplacian" so that K = L + N, where L is the "usual"
  %     cotangent Laplacian and N computes normal derivatives at boundary
  %     vertices.
  %

  if size(Ele,2) == 3
    assert(isempty(find_ears(Ele)), ...
      'Mesh should not have ears, preprocess Ele with flip_ears');
  end

  % http://www.cs.toronto.edu/~jacobson/images/error-in-linear-subspace-design-for-real-time-shape-deformation-2017-wang-et-al.pdf
  use_paper_version = false;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'PaperVersion'}, ...
    {'use_paper_version'});
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


  if use_paper_version
    L = cotmatrix(V,Ele);
    [DD,TE] = normal_derivative(V,Ele);
    AT = sparse( ...
      TE(:), ...
      repmat(1:size(TE,1),1,size(TE,2))', ...
      1,size(V,1),size(TE,1));
    [~,C] = on_boundary(Ele);
    DD(~C,:) = 0;
    N = (0.5*AT*DD);
    K = L + N;
    M = massmatrix(V,Ele);
    A = K' * (M\K);
  else
    [Mcr,E,EMAP] = crouzeix_raviart_massmatrix(V,Ele);
    [Lcr] = crouzeix_raviart_cotmatrix(V,Ele);
    [~,C] = on_boundary(Ele);
    % Ad  #E by #V Edge-vertex incidence matrix
    Ad = sparse(E(:),repmat(1:size(E,1),1,size(E,2))',1,size(V,1),size(E,1))';
    De = diag(sparse(sum(Ad,2)));
    % Invert mass matrix
    iMcr = diag(sparse(1./diag(Mcr)));
    % kill boundary edges
    iMcr(EMAP(C),EMAP(C)) = 0;
    Le = Lcr*(De\Ad);
    A = Le'*(iMcr*Le);
  end

  if nargin<3 || isempty(b)
    W = [];
  else
    bc = eye(numel(b));
    W = min_quad_with_fixed(A,[],b,bc);
  end
end
