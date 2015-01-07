function [b,K] = arap_rhs(varargin)
  % ARAP_RHS build right-hand side of global poisson solve for various ARAP
  % energies. For examples see Sorkine and Alexa's ARAP (see equation 8 and 9):
  %     ∑ wij * 0.5 * (Ri + Rj) * (Vi - Vj)
  %   j∈N(i)
  % 
  % b = arap_rhs(V,F,R)
  % [b,K] = arap_rhs(V,F,R,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by dim list of initial domain positions
  %   F  #F by 3 list of triangle indices into V
  %   R  dim by dim by #V list of rotations, if R is empty then b is set to []
  %   Optional:
  %     'Energy'
  %       followed by a string specifying which arap energy definition to use.
  %       One of the following:
  %         'spokes'  "As-rigid-as-possible Surface Modeling" by [Sorkine and
  %           Alexa 2007], rotations defined at vertices affecting incident
  %           edges, default
  %         'elements'  "A local-global approach to mesh parameterization" by
  %           [Liu et al.  2010] or "A simple geometric model for elastic
  %           deformation" by [Chao et al.  2010], rotations defined at
  %           elements (triangles or tets) 
  %         'spokes-and-rims'  Adapted version of "As-rigid-as-possible Surface
  %           Modeling" by [Sorkine and Alexa 2007] presented in section 4.2 of
  %           or "A simple geometric model for elastic deformation" by [Chao et
  %           al.  2010], rotations defined at vertices affecting incident
  %           edges and opposite edges
  %     and only K is computed
  %      'Dim'  solve dimension {size(V,2)}
  % Output:
  %   b  #V by dim right hand side
  %   K  #V*dim by #(F|V)*dim*dim matrix such that: 
  %     b = K * reshape(permute(R,[3 1 2]),size(V,1)*size(V,2)*size(V,2),1);
  %   

  % default is Sorkine and Alexa style local rigidity energy
  energy = 'spokes';
  V = varargin{1};
  F = varargin{2};
  R = varargin{3};
  % number of vertices
  n = size(V,1);
  % number of elements
  m = size(F,1);
  % simplex size
  simplex_size = size(F,2);
  assert(simplex_size == 3 || simplex_size == 4);
  % number of dimensions
  Vdim = size(V,2);
  dim = Vdim;

  ii = 4;
  while(ii <= nargin)
    switch varargin{ii}
    case 'Energy'
      ii = ii + 1;
      assert(ii<=nargin);
      energy = varargin{ii};
    case 'Dim'
      ii = ii + 1;
      assert(ii<=nargin);
      dim = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii + 1;
  end

  % number of rotations
  switch energy
  case 'spokes'
    nr = size(V,1);
  case 'spokes-and-rims'
    nr = size(V,1);
  case 'elements'
    nr = size(F,1);
  end
  
  KX = arap_linear_block(V,F,1,'Energy',energy);
  KY = arap_linear_block(V,F,2,'Energy',energy);
  Z = sparse(size(V,1),nr);
  if Vdim == 2
    K = [ ...
      KX Z KY Z; ...
      Z KX Z KY];
    assert(dim == 2);
  elseif Vdim == 3
    KZ = arap_linear_block(V,F,3,'Energy',energy);
    K = [ ...
             KX Z Z  KY Z  Z  KZ Z  Z ; ...
             Z KX Z  Z  KY Z  Z  KZ Z];
    if dim == 3
      K = [K;Z Z  KX Z  Z  KY Z  Z  KZ ];
    else
      assert(dim == 2);
    end
  end

  if(~isempty(R))
    assert(dim == size(R,1));
    assert(dim == size(R,2));
    assert(nr == size(R,3));
    % collect rotations into a single column
    Rcol = reshape(permute(R,[3 1 2]),nr*dim*dim,1);
    b = K * Rcol;
    b = reshape(b,[size(V,1) dim]);
  else
    b = [];
  end

  % Notes on construction for 'spokes' energy
  %
  % Bi =   ∑ wij * 0.5 * (Ri + Rj) * (Vi - Vj)
  %      j∈N(i)
  %
  % where Bi, Vi, and Vj are dim-length vectors, and Ri and Rj are
  % dim by dim rotation matrices
  %
  % Bi' =   ∑ wij * 0.5 * (Vi'-Vj') * (Ri' + Rj')
  %       j∈N(i)
  %
  % Bi(x) =   ∑ wij * 0.5 * ∑ (Vi(y)-Vj(y)) * (Ri(x,y) + Rj(x,y))
  %          j∈N(i)        y∈dim
  %
  % Where Bi(x) is getting the xth coordinate of Bi and Ri(x,y) is getting the
  % entry of Ri at row x and col y
  % 

end
