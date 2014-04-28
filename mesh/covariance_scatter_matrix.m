function [CSM] = covariance_scatter_matrix(varargin)
  % COVARIANCE_SCATTER_MATRIX builds a scatter matrix for constructing the
  % covariance matrices used in ARAP energy minimization. 
  %
  % CSM = covariance_scatter_matrix(V,F)
  % CSM = covariance_scatter_matrix(V,F,'ParameterName',ParameterValue)
  % 
  % Inputs:
  %   V  #V by dim list of initial domain positions
  %   F  #F by #simplex size list of triangle indices into V
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
  %   CSM dim*nr by dim*n sparse matrix containing special laplacians along the
  %     diagonal so that when multiplied by repmat(U,dim,1) gives covariance 
  %     matrix elements, can be used to speed up covariance matrix computation.
  %     Where nr is the number of best fit rotations which depends on the
  %     choice of energy
  % Example:
  %  % compute covariance matrices of a mesh
  %  n = size(V,1);
  %  dim = size(V,2);
  %  CSM = covariance_scatter_matrix(V,F,'Energy','spokes')
  %  S = CSM*repmat(U,dim,1);
  %  % dim by dim by n list of covariance matrices
  %  S = permute(reshape(S,[n dim dim]),[2 3 1]);
  %

  % default is Sorkine and Alexa style local rigidity energy
  energy = 'spokes';
  V = varargin{1};
  F = varargin{2};
  % number of vertices
  n = size(V,1);
  % number of elements
  m = size(F,1);
  % number of dimensions
  dim = size(V,2);

  ii = 3;
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
  if dim == 2
    CSM = [KX Z;Z KY]';
  elseif dim == 3
    KZ = arap_linear_block(V,F,3,'Energy',energy);
    CSM = [KX Z Z;Z KY Z;Z Z KZ]';
  end

end
