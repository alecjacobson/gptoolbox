function [CV,CF,CJ,CI] = edge_cylinders(PV,PE,varargin);
  % EDGE_CYLINDERS Generate a cylinder mesh at each edge, no accounting for
  % overlaps.
  %
  % [CV,CF,CJ,CI] = edge_cylinders(P,PE,varargin);
  %
  % Inputs:
  %   PV  #PV by 3 list of 3d edge vertex positions
  %   PE  #PE by 2 list of edge indices into PV
  %   Optional:
  %   'Thickness' followed by diameter thickness of wire {0.1 average edge length}
  %   'PolySize'  followed by number of sides on each wire (e.g., 4 would produce
  % Outputs:
  %   CV  #CV by 3 list of cylinder vertices
  %   CF  #CF by 3 list of cylinder face indices into CV
  %   CJ  #CF list of indices into PE
  %   CI  #CV list of indices into PV
  %
  % See also: wire_mesh
  %

  function mod = quatmod( q )
    %  QUATMOD Calculate the modulus of a quaternion.
    %   N = QUATMOD( Q ) calculates the modulus, N, for a given quaternion, Q.  
    %   Input Q is an M-by-4 matrix containing M quaternions.  N returns a 
    %   column vector of M moduli.  Each element of Q must be a real number.  
    %   Additionally, Q has its scalar number as the first column.
    %
    %   Examples:
    %
    %   Determine the modulus of q = [1 0 0 0]:
    %      mod = quatmod([1 0 0 0])
    %
    %   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMULTIPLY, QUATNORM,
    %   QUATNORMALIZE, QUATROTATE.
    
    %   Copyright 2000-2011 The MathWorks, Inc.
    
    if any(~isreal(q(:)))
        error(message('aero:quatnorm:isNotReal'));
    end
    
    if (size(q,2) ~= 4)
        error(message('aero:quatnorm:wrongDimension'));
    end
    
    for index = size(q,1):-1:1
        mod(index,:) = norm(q(index,:),2);
    end
  end

  function qout = quatnormalize( q )
    %  QUATNORMALIZE Normalize a quaternion.
    %   N = QUATNORMALIZE( Q ) calculates the normalized quaternion, N, for a
    %   given quaternion, Q.  Input Q is an M-by-4 matrix containing M
    %   quaternions.  N returns an M-by-4 matrix of normalized quaternions.
    %   Each element of Q must be a real number.  Additionally, Q has its
    %   scalar number as the first column.
    %
    %   Examples:
    %
    %   Normalize q = [1 0 1 0]:
    %      normal = quatnormalize([1 0 1 0])
    %
    %   See also QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD, QUATMULTIPLY, 
    %   QUATNORM, QUATROTATE.
    
    %   Copyright 2000-2005 The MathWorks, Inc.
    
    qout = q./(quatmod( q )* ones(1,4));
  end
  function dcm = quat2dcm( q )
    %  QUAT2DCM Convert quaternion to direction cosine matrix.
    %   N = QUAT2DCM( Q ) calculates the direction cosine matrix, N, for a
    %   given quaternion, Q.  Input Q is an M-by-4 matrix containing M
    %   quaternions.  N returns a 3-by-3-by-M matrix of direction cosine
    %   matrices.  The direction cosine matrix performs the coordinate
    %   transformation of a vector in inertial axes to a vector in body axes.
    %   Each element of Q must be a real number.  Additionally, Q has its
    %   scalar number as the first column.
    %
    %   Examples:
    %
    %   Determine the direction cosine matrix from q = [1 0 1 0]:
    %      dcm = quat2dcm([1 0 1 0])
    %
    %   Determine the direction cosine matrices from multiple quaternions:
    %      q = [1 0 1 0; 1 0.5 0.3 0.1];
    %      dcm = quat2dcm(q)
    %
    %   See also ANGLE2DCM, DCM2ANGLE, DCM2QUAT, ANGLE2QUAT, QUAT2ANGLE, QUATROTATE.
    
    %   Copyright 2000-2007 The MathWorks, Inc.
    
    qin = quatnormalize( q );
    
    dcm = zeros(3,3,size(qin,1));
    
    dcm(1,1,:) = qin(:,1).^2 + qin(:,2).^2 - qin(:,3).^2 - qin(:,4).^2;
    dcm(1,2,:) = 2.*(qin(:,2).*qin(:,3) + qin(:,1).*qin(:,4));
    dcm(1,3,:) = 2.*(qin(:,2).*qin(:,4) - qin(:,1).*qin(:,3));
    dcm(2,1,:) = 2.*(qin(:,2).*qin(:,3) - qin(:,1).*qin(:,4));
    dcm(2,2,:) = qin(:,1).^2 - qin(:,2).^2 + qin(:,3).^2 - qin(:,4).^2;
    dcm(2,3,:) = 2.*(qin(:,3).*qin(:,4) + qin(:,1).*qin(:,2));
    dcm(3,1,:) = 2.*(qin(:,2).*qin(:,4) + qin(:,1).*qin(:,3));
    dcm(3,2,:) = 2.*(qin(:,3).*qin(:,4) - qin(:,1).*qin(:,2));
    dcm(3,3,:) = qin(:,1).^2 - qin(:,2).^2 - qin(:,3).^2 + qin(:,4).^2;
    end

  function qout = quatrotate( q, r )
    %  QUATROTATE Rotate a vector by a quaternion.
    %   N = QUATROTATE( Q, R ) calculates the rotated vector, N, for a
    %   quaternion, Q, and a vector, R.  Q is either a M-by-4 matrix
    %   containing M quaternions or a single 1-by4 quaternion.  R
    %   is either a M-by-3 matrix or a single 1-by-3 vector.  N returns an
    %   M-by-3 matrix of rotated vectors.  Each element of Q and R must be a
    %   real number.  Additionally, Q has its scalar number as the first column.
    %
    %   Examples:
    %
    %      q = [1 0 1 0];
    %      r = [1 1 1];
    %      n = quatrotate( q, r ) 
    %
    %      q = [1 0 1 0; 1 0.5 0.3 0.1];
    %      r = [1 1 1];
    %      n = quatrotate( q, r ) 
    %
    %      q = [1 0 1 0];
    %      r = [1 1 1; 2 3 4];
    %      n = quatrotate( q, r ) 
    %
    %      q = [1 0 1 0; 1 0.5 0.3 0.1];
    %      r = [1 1 1; 2 3 4];
    %      n = quatrotate( q, r ) 
    %
    %   See also QUAT2DCM, QUATCONJ, QUATDIVIDE, QUATINV, QUATMOD,
    %   QUATMULTIPLY, QUATNORM, QUATNORMALIZE.

    %   Copyright 2000-2010 The MathWorks, Inc.
        
    if any(~isreal(q(:)))
        error(message('aero:quatrotate:isNotReal1'));
    end

    if (size(q,2) ~= 4)
        error(message('aero:quatrotate:wrongDimension1'));
    end
        
    if any(~isreal(r(:)))
        error(message('aero:quatrotate:isNotReal2'));
    end

    if (size(r,2) ~= 3)
        error(message('aero:quatrotate:wrongDimension2'));
    end

    if (size(r,1) ~= size(q,1) && ~( size(r,1) == 1 || size(q,1) == 1))
        error(message('aero:quatrotate:wrongDimension3'));
    end

    dcm = quat2dcm(q);

    if ( size(q,1) == 1 ) 
        % Q is 1-by-4
        qout = (dcm*r')';
    elseif (size(r,1) == 1) 
        % R is 1-by-3
        for i = size(q,1):-1:1
            qout(i,:) = (dcm(:,:,i)*r')';
        end
    else
        % Q is M-by-4 and R is M-by-3
        for i = size(q,1):-1:1
            qout(i,:) = (dcm(:,:,i)*r(i,:)')';
        end
    end
  end


  thickness = 0.1*mean(normrow(PV(PE(:,1),:)-PV(PE(:,2))));
  poly = 4;

  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'PolySize','Thickness'}, ...
    {'poly','thickness'});
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

  % All unique undirected edges
  ne = size(PE,1);
  % edge lengths
  Evec = PV(PE(:,2),:)-PV(PE(:,1),:);
  l = sqrt(sum(Evec.^2,2));
  [w,a] = axisanglebetween(Evec,repmat([0 0 1],size(Evec,1),1));
  bad = any(isnan(w),2);
  w(bad,1) = 1; w(bad,2) = 0; w(bad,3) = 0;
  a(bad) = (Evec(bad,3)<0)*pi;

  % twice as thick
  [x,y,z] = cylinder(thickness/2,poly);
  [CF,CV] = surf2patch(x,y,z,'triangles');
  CI = repmat([1;2],poly,1);
  [CV,I,J] = remove_duplicate_vertices(CV,eps);
  CI = CI(I);
  CF = J(CF);
  CF = [CF;fill_holes(CV,CF)];
  ncv = size(CV,1);
  ne = size(PE,1);
  CJ = reshape(repmat(1:ne,size(CF,1),1),ne*size(CF,1),1);
  CF = reshape(permute(bsxfun(@plus,permute((0:ne-1)*size(CV,1),[1 3 2]),CF),[1 3 2]),ne*size(CF,1),3);
  CV = reshape(permute( bsxfun(@times,permute([ones(ne,2) l],[3 2 1]),CV),[1 3 2]),ne*ncv,3);
  CI = PE(sub2ind(size(PE),repmat(1:ne,size(CI,1),1),repmat(CI,1,ne)));
  CI = CI(:);
  R  = reshape(permute(repmat(permute(axisangle2quat(w,a),[3 2 1]),ncv,1),[1 3 2]),ne*ncv,4);
  T  = reshape(permute(        repmat(permute(PV(PE(:,1),:),[3 2 1]),ncv,1),[1 3 2]),ne*ncv,3);
  CV = quatrotate(R,CV)+T;

end
