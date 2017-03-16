function [Kd] = arap_linear_block(varargin)
  % ARAP_LINEAR_BLOCK constructs a block of the matrix which constructs the
  % linear terms of a given arap energy. When treating rotations as knowns
  % (arranged in a column) then this constructs Kd of K such that the linear
  % portion of the energy is as a column:
  %   K * R = [Kx Z  ... Ky Z  ... 
  %            Z  Kx ... Z  Ky ... 
  %            ... ]
  % These blocks are also used to build the "covariance scatter matrices". Here
  % we want to build a scatter matrix that multiplies against positions
  % (treated as known) producing covariance matrices to fit each rotation.
  % Notice that in the case of the RHS of the poisson solve the rotations are
  % known and the positions unknown, and vice versa for rotation fitting. These
  % linear block just relate the rotations to the positions, linearly in each.
  %
  % Kd = arap_linear_block(V,F,d)
  % Kd = arap_linear_block(V,F,d,'ParameterName','ParameterValue,...)
  % 
  % Inputs:
  %   V  #V by dim list of initial domain positions
  %   F  #F by #simplex size list of triangle indices into V
  %   d  coordinate of linear constructor to build
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
  % Outputs:
  %   Kd  #V by #V/#F block of the linear constructor matrix corresponding to
  %     coordinate d
  %
  % See also: arap, arap_rhs, covariance_scatter_matrix
  %

  % default is Sorkine and Alexa style local rigidity energy
  energy = 'spokes';
  V = varargin{1};
  F = varargin{2};
  d = varargin{3};
  % number of vertices
  n = size(V,1);
  % number of elements
  m = size(F,1);
  % simplex size
  simplex_size = size(F,2);
  assert(simplex_size == 3 || simplex_size == 4);
  % number of dimensions
  dim = size(V,2);
  assert(d <= dim);

  ii = 4;
  while(ii <= nargin)
    switch varargin{ii}
    case 'Energy'
      ii = ii + 1;
      assert(ii<=nargin);
      energy = varargin{ii};
    otherwise
      error(['Unsupported parameter: ' varargin{ii}]);
    end
    ii = ii + 1;
  end

  switch energy
  case 'spokes'
    Kd = spokes_linear_block(V,F,d);
  case 'spokes-and-rims'
    Kd = spokes_and_rims_linear_block(V,F,d);
  case 'elements'
    Kd = elements_linear_block(V,F,d);
  otherwise
    error(['Unsupported energy type: ' energy]);
  end

  function K = spokes_linear_block(V,F,d)
    % Computes a matrix K such that V'* K * R computes
    %  ∑ wij * 0.5 * (V(i,d)-V(j,d)) * (Ri + Rj)
    % j∈N(i)
    % 
    % Inputs:
    %   V  #V by dim list of coordinates
    %   F  #F by 3 list of triangle indices into V
    %   d  index into columns of V
    % Output:
    %   K  #V by #V matrix
    %
    E = edges(F);
    % Build upper part of adjacency matrix where instead of a 1 for edge from i
    % to j we have the difference of position in dimension d
    A = sparse(E(:,1),E(:,2),V(E(:,1),d)-V(E(:,2),d),size(V,1),size(V,1));
    % anti-symmetric, or considers direction of edges
    A = A-A';
    % Multiply with cotangent weights (don't worry about diagonal begin wrong
    % since in A it's all zeros
    L = cotmatrix(V,F);
    K = L.*A;
    % correct the diagonal (notice that the sign is positive
    K = K + diag(sum(K,2));
    K = 0.5*K;
  end

  function K = spokes_and_rims_linear_block(V,F,d)
    % Computes a matrix K such that V' * K * R computes
    %  ∑    -2*(cot(aij) + cot(bij) * (V(i,d)-V(j,d)) * (Ri + Rj) +
    %       -2*cot(aij) * (V(i,d)-V(j,d)) * Raij +
    %       -2*cot(bij) * (V(i,d)-V(j,d)) * Rbij
    % j∈N(i)
    % 
    % where:          vj
    %              /  |  \
    %             /   |   \
    %            /    |    \
    %           /     |     \
    %          aij    |    bij
    %           \     |     /
    %            \    |    /
    %             \fij|gij/
    %              \  |  /
    %                 vi
    % 
    % Inputs:
    %   V  #V by dim list of coordinates
    %   F  #F by 3 list of triangle indices into V
    %   d  index into columns of V
    % Output:
    %   K  #V by #F matrix
    %
    if simplex_size == 3
      % triangles
      C = cotangent(V,F);
      i1 = F(:,1); i2 = F(:,2); i3 = F(:,3);
      I = [i1;i2;i2;i3;i3;i1;i1;i2;i3];
      J = [i2;i1;i3;i2;i1;i3;i1;i2;i3];
      v = [ ...
         C(:,3).*(V(i1,d)-V(i2,d)) + C(:,2).*(V(i1,d)-V(i3,d)); ... 
        -C(:,3).*(V(i1,d)-V(i2,d)) + C(:,1).*(V(i2,d)-V(i3,d)); ... 
         C(:,1).*(V(i2,d)-V(i3,d)) + C(:,3).*(V(i2,d)-V(i1,d)); ... 
        -C(:,1).*(V(i2,d)-V(i3,d)) + C(:,2).*(V(i3,d)-V(i1,d)); ... 
         C(:,2).*(V(i3,d)-V(i1,d)) + C(:,1).*(V(i3,d)-V(i2,d)); ... 
        -C(:,2).*(V(i3,d)-V(i1,d)) + C(:,3).*(V(i1,d)-V(i2,d)); ... 
         ... % diagonal
        C(:,3).*(V(i1,d)-V(i2,d)) - C(:,2).*(V(i3,d)-V(i1,d)); ... 
        C(:,1).*(V(i2,d)-V(i3,d)) - C(:,3).*(V(i1,d)-V(i2,d)); ... 
        C(:,2).*(V(i3,d)-V(i1,d)) - C(:,1).*(V(i2,d)-V(i3,d)); ... 
        ];
      % construct and divide by 3 so laplacian can be used as is
      K = sparse(I,J,v,n,n)/3;
    elseif simplex_size == 4
      % tetrahedra
      assert(false)
    end
  end

  function K = elements_linear_block(V,F,d)
    % Computes a matrix K such that V' * K * R computes
    %  ∑    -2*cot(aij) * (V(i,d)-V(j,d)) * Raij +
    %       -2*cot(bij) * (V(i,d)-V(j,d)) * Rbij
    % j∈N(i)
    % 
    % where:        vj
    %              / |\
    %             /  | \
    %            aij | bij
    %             \  | /
    %              \ |/
    %               vi
    % 
    % Inputs:
    %   V  #V by dim list of coordinates
    %   F  #F by 3 list of triangle indices into V
    %   d  index into columns of V
    % Output:
    %   K  #V by #F matrix
    %
    if simplex_size == 3
      % triangles, #T by 3
      C = cotangent(V,F);
      i1 = F(:,1); i2 = F(:,2); i3 = F(:,3);
      I = [i2;i3;i3;i1;i1;i2];
      J = repmat((1:m)',3*2,1);
      v = [ ...
         C(:,1).*(V(i2,d)-V(i3,d)); ...
        -C(:,1).*(V(i2,d)-V(i3,d)); ...
         C(:,2).*(V(i3,d)-V(i1,d)); ...
        -C(:,2).*(V(i3,d)-V(i1,d)); ...
         C(:,3).*(V(i1,d)-V(i2,d)); ...
        -C(:,3).*(V(i1,d)-V(i2,d))];
      K = sparse(I,J,v,n,m);
    elseif simplex_size == 4
      % tetrahedra, #T by 6
      C = cotangent(V,F);
      i1 = F(:,1); i2 = F(:,2); i3 = F(:,3); i4 = F(:,4);
      I = [i2;i3;i3;i1;i1;i2;i4;i1;i4;i2;i4;i3];
      J = repmat((1:m)',6*2,1);
      v = [ ...
         C(:,1).*(V(i2,d)-V(i3,d)); ...
        -C(:,1).*(V(i2,d)-V(i3,d)); ...
         C(:,2).*(V(i3,d)-V(i1,d)); ...
        -C(:,2).*(V(i3,d)-V(i1,d)); ...
         C(:,3).*(V(i1,d)-V(i2,d)); ...
        -C(:,3).*(V(i1,d)-V(i2,d)); ...
         C(:,4).*(V(i4,d)-V(i1,d)); ...
        -C(:,4).*(V(i4,d)-V(i1,d)); ...
         C(:,5).*(V(i4,d)-V(i2,d)); ...
        -C(:,5).*(V(i4,d)-V(i2,d)); ...
         C(:,6).*(V(i4,d)-V(i3,d)); ...
        -C(:,6).*(V(i4,d)-V(i3,d))];
      K = sparse(I,J,v,n,m);
    end
  end

end
