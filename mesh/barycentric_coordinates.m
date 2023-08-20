function B = barycentric_coordinates(P,varargin)
  % BARYCENTRIC_COORDINATES Computes barycentric coordinates of point p in
  % simplex (v1,v2,v3).
  % 
  % B = barycentric_coordinates(P,V1,V2,V3, ...)
  %
  % Inputs:
  %   P  #P by dim list of query point locations
  %   V1  #P by dim list of simplex corner locations
  %   V2  #P by dim list of simplex corner locations
  %   ...
  %   Optional:
  %     'Project'  followed by whether to project P onto the hyperplane spanned
  %     by V1,V2,... {false}
  % Outputs:
  %   B  #P by dim+1 list of barycentric coordinates
  %

  ss = 0;
  for varg = varargin
    if ischar(varg{1})
      break;
    end
    assert(size(P,1) == size(varg{1},1), 'All inputs should be same length');
    assert(size(P,2) == size(varg{1},2), 'All inputs should be same dimension');
    ss = ss+1;
  end

  project = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Project'}, ...
    {'project'});
  v = ss+1;
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
  dim = size(P,2);
  project = project || dim+1==ss;

  function v = volume(ad,bd,cd)
    r =[bd(:,2).*cd(:,3)-bd(:,3).*cd(:,2), ...
        bd(:,3).*cd(:,1)-bd(:,1).*cd(:,3), ...
        bd(:,1).*cd(:,2)-bd(:,2).*cd(:,1)];
    v = -sum(ad.*r,2)./6;
  end

  n = size(P,1);
  if project && ss<4
    switch ss
    case 4
      assert(dim == 3,'Project=true only implemented for triangles');
      V1 = varargin{1};
      V2 = varargin{2};
      V3 = varargin{3};
      V4 = varargin{4};
      V1P = V1-P;
      V2P = V2-P;
      V3P = V3-P;
      V4P = V4-P;
      %T = bsxfun(@plus,(1:n)',(0:3)*n);
      A1 = volume(V2P,V4P,V3P);
      A2 = volume(V1P,V3P,V4P);
      A3 = volume(V1P,V4P,V2P);
      A4 = volume(V1P,V2P,V3P);
      A  = volume(V1-V4,V2-V4,V3-V4);
      if dim>3 && max(abs(sum([A1 A2 A3 A4],2)-A))>1e-14
        warning('Possibly negative coordinates. Not supported in dim~=3');
      end
      B = bsxfun(@rdivide,[A1 A2 A3 A4],A);
    case 3
      V1 = varargin{1};
      V2 = varargin{2};
      V3 = varargin{3};
      switch dim
      case 2
        A1 = doublearea([ P;V2;V3],[1:n;n+[1:n;n+(1:n)]]');
        A2 = doublearea([V1; P;V3],[1:n;n+[1:n;n+(1:n)]]');
        A3 = doublearea([V1;V2; P],[1:n;n+[1:n;n+(1:n)]]');
        A  = doublearea([V1;V2;V3],[1:n;n+[1:n;n+(1:n)]]');
        B = bsxfun(@rdivide,[A1 A2 A3],A);
      case 3
        u = V2-V1;
        v = V3-V1;
        N = cross(u,v,2);
        w = P-V1;
        N2 = sum(N.^2,2);
        gamma = sum(cross(u,w,2).*N,2)./N2;
        beta = sum(cross(w,v,2).*N,2)./N2;
        alpha = 1-gamma-beta;
        B = [alpha beta gamma];
      otherwise
        error('Cannot use project=true for dim>2');
      end
    case 2
      V1 = varargin{1};
      V2 = varargin{2};
      A1 = edge_lengths([P;V2],[1:n;n+[1:n]]');
      A2 = edge_lengths([V1;P],[1:n;n+[1:n]]');
      A = edge_lengths([V1;V2],[1:n;n+[1:n]]');
      if dim>1 && max(abs(sum([A1 A2],2)-A))>1e-14
        warning('Possibly negative coordinates. Not supported in dim~=1');
      end
      B = bsxfun(@rdivide,[A1 A2],A);
    otherwise
      error(sprintf('%d-simplices not supported',ss));
    end
  else
    B = zeros(size(P,1),ss);
    Ai = [zeros(ss,dim) ones(ss,1)];
    for i = 1:size(P,1)
      for s = 1:ss
        Ai(s,1:dim) = varargin{s}(i,:);
      end
      B(i,:) = [P(i,:) 1] / Ai;
    end
  end
end
