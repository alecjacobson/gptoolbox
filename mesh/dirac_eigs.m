function [eVec,eVal] = dirac_eigs(V,F,k,varargin)
  % DIRAC_EIGS Construct the eigenvalue and eigenvectors of the discrete
  % quaternionic dirac operator for a 3d triangle mesh (as described by "Spin
  % Transformations of Discrete Surfaces" [Crane et al. 2011]).
  % 
  % [eVec,eVal] = dirac_eigs(V,F,k,varargin)
  %
  % Inputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   k  number of eigenmodes to compute
  %   Optional:
  %     'Relative' followed by wehter to use the "relative dirac operator"
  %        (as described by "A Dirac Operator for Extrinsic Shape Analysis"
  %        [Liu, Jacobson, and Crane. 2017]).) {false}.
  %     'Tao'  followed by how much to blend between "relative dirac operator"
  %     and Laplace opeator. Only used if {'Relative',true}. {0.999999}
  % Outputs:
  %   eVec  #V by 4 by k list of eigenvectors. Each 4x1 subblock
  %     represents a quaternion (w,x,y,z)
  %   eVal  k by 1 list of eigenvalues
  %

  function [eVal, eVec] = QuaternionEigs(mat, massMat, numEigs)
    opts.issym = 1;
    opts.isreal = 1;
    opt.tol = 1e-7;
    sigma = 1e-6;
    [eVec, eVal] = eigs(mat + sigma.*speye(size(mat,1),size(mat,2)),...
      massMat, numEigs*4, 'sm', opt);
    [~ ,I] = sort(diag(abs(eVal))); % sort
    eVal = sum(eVal,2);
    eVal = eVal(I);
    eVal = real(eVal) - sigma;
    eVec = real(eVec(:,I));
  end

  % default values
  relative = false;
  tao = 0.999999;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Relative','Tao'}, ...
    {'relative','tao'});
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
    
  if relative
    D = relative_dirac_operator(V,F);
  else
    D = dirac_operator(V,F);
  end

  % 1. Compute Relative Dirac operator
  if relative
    D = relative_dirac_operator(V,F);
  else
    D = dirac_operator(V,F);
  end
  
  % 2. Solve the eigenvalue problem
  A = doublearea(V,F)/2;
  MF= kroneye(diag(sparse(A)),4);
  Dsquare = D'*MF*D; % the square of the operator
  % average from vertices to faces
  B = kroneye(sparse(repmat(1:size(F,1),3,1)',F,1/3,size(F,1),size(V,1)),4);
  massMat = B'*MF*D; % "mass matrix"
  massMat = (massMat'+massMat)/2; % Symmetrize mass matrix

  if relative && tao<1
    L = kroneye(cotmatrix(V,F),4);
    M = kroneye(massmatrix(V,F),4);
    Dsquare = Dsquare + tao*(L-Dsquare);
    massMat = massMat + tao*(M-massMat);

  end

  [eVal, eVec] = QuaternionEigs(Dsquare, massMat, k); 

  % 3. normalize eigenvectors
  vertArea = full(diag(massmatrix(V,F,'barycentric')));
  vertArea4D = zeros(4*size(V,1),1);
  vertArea4D(1:4:end) = vertArea;
  vertArea4D(2:4:end) = vertArea;
  vertArea4D(3:4:end) = vertArea;
  vertArea4D(4:4:end) = vertArea;
  ALambda = repmat(vertArea4D, [1, size(eVec,2)]) .* eVec;
  Z = zeros(4, size(eVec,2));
  Z(1,:) = sum(ALambda(1:4:end,:));
  Z(2,:) = sum(ALambda(2:4:end,:));
  Z(3,:) = sum(ALambda(3:4:end,:));
  Z(4,:) = sum(ALambda(4:4:end,:));

  % Calculate quaternion inverse
  mag = sum(Z.^2)';
  Z(2,:) = -Z(2,:); Z(3,:) = -Z(3,:); Z(4,:) = -Z(4,:);
  Z = Z' ./ repmat(mag, 1, 4);

  % reshape so that modes are in the 3rd array dimension
  Z = permute(Z,[3 2 1]);
  n = size(V,1);
  eVec = permute(reshape(eVec,[4 n 4*k]),[2 1 3]);

  % % Vectorized quatmultiply
  % q = eVec;
  % r = Z;
  % vec = [q(:,1,:).*repmat(r(:,2,:),size(q,1),1) q(:,1,:).*repmat(r(:,3,:),size(q,1),1) q(:,1,:).*repmat(r(:,4,:),size(q,1),1)] + ...
  %       [repmat(r(:,1,:),size(q,1),1).*q(:,2,:) repmat(r(:,1,:),size(q,1),1).*q(:,3,:) repmat(r(:,1,:),size(q,1),1).*q(:,4,:)]+...
  %       [q(:,3,:).*repmat(r(:,4,:),size(q,1),1)-q(:,4,:).*repmat(r(:,3,:),size(q,1),1) ...
  %         q(:,4,:).*repmat(r(:,2,:),size(q,1),1)-q(:,2,:).*repmat(r(:,4,:),size(q,1),1) ...
  %         q(:,2,:).*repmat(r(:,3,:),size(q,1),1)-q(:,3,:).*repmat(r(:,2,:),size(q,1),1)];
  % scalar = q(:,1,:).*repmat(r(:,1,:),size(q,1),1) - q(:,2,:).*repmat(r(:,2,:),size(q,1),1) - q(:,3,:).*repmat(r(:,3,:),size(q,1),1) - q(:,4,:).*repmat(r(:,4,:),size(q,1),1);
  % normEVec = cat(2,scalar,vec);

  % Vectorized quatmultiply
  q = eVec;
  r = Z;
  vec = [q(:,1,:).*r(:,2,:) q(:,1,:).*r(:,3,:) q(:,1,:).*r(:,4,:)] + ...
        [r(:,1,:).*q(:,2,:) r(:,1,:).*q(:,3,:) r(:,1,:).*q(:,4,:)]+...
        [q(:,3,:).*r(:,4,:)-q(:,4,:).*r(:,3,:) ...
         q(:,4,:).*r(:,2,:)-q(:,2,:).*r(:,4,:) ...
         q(:,2,:).*r(:,3,:)-q(:,3,:).*r(:,2,:)];
  scalar = q(:,1,:).*r(:,1,:) - q(:,2,:).*r(:,2,:) - q(:,3,:).*r(:,3,:) - q(:,4,:).*r(:,4,:);
  normEVec = cat(2,scalar,vec);

  % reshape so that modes are back in the 2nd array dimension
  normEVec = reshape(permute(normEVec,[2 1 3]),[4*n 4*k]);
  normEVec = normEVec(:,1:4:end);
  normEVec = normalizerow(normEVec')';

  % 4. Output
  eVal = eVal(1:4:end);
  eVec = normEVec;
  eVec = permute(cat(3, ...
  eVec(1:4:end,:),eVec(2:4:end,:),eVec(3:4:end,:),eVec(4:4:end,:)),[1 3 2]);

end
