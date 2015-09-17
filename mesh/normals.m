function [ N ] = normals(V,F,varargin)
  % NORMALS Compute *unnormalized* normals per face
  %
  % N = normals(V,F)
  % N = normals(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %  V  #V x 3 matrix of vertex coordinates
  %  F  #F x 3  matrix of indices of triangle corners
  %  Optional:
  %    'Stable' followed by whether to compute normals in a way stable with
  %      respect to vertex order: constant factor more expensive {false}
  %    'UseSVD' followed by whether to use SVD, slow {false}
  % Output:
  %  N  #F x 3 list of face normals
  %

  function D = sum3(A,B,C)
    % SUM3 Entrywise sum of three matrices in a stable way: sorting entries by
    % value and summing 
    shape = size(A);
    ABC = [A(:) B(:) C(:)];
    [~,I] = sort(abs(ABC),2,'descend');
    sABC = reshape( ...
      ABC(sub2ind(size(ABC),repmat(1:size(ABC,1),size(ABC,2),1)',I)),size(ABC));
    D = reshape(sum(sABC,2),shape);
  end

  stable = false;
  use_svd = false;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Stable','UseSVD'},{'stable','use_svd'});
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
  
  p1 = V(F(:,1),:);
  p2 = V(F(:,2),:);
  p3 = V(F(:,3),:);
  
  if use_svd
    N = zeros(size(F,1),3);
    BC = barycenter(V,F);
    for f = 1:size(F,1)
      Uf = bsxfun(@minus,V(F(f,:),:),BC(f,:));
      [~,~,sV] = svd(Uf);
      N(f,:) = sV(:,3);
    end
    NN = normals(V,F,'UseSVD',false);
    N(sum(N.*NN,2)<0,:) = N(sum(N.*NN,2)<0,:)*-1;
  else
    % ,2 is necessary because this will produce the wrong result if there are
    % exactly 3 faces
    N1 = cross(p2 - p1, p3 - p1,2);
    if stable
      N2 = cross(p3 - p2, p1 - p2,2);
      N3 = cross(p1 - p3, p2 - p3,2);
      N = sum3(N1,N2,N3)/3;
    else
      N = N1;
    end
  end


end

