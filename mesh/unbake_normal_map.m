function [V,F,TV,TF,obj] = unbake_normal_map(V,F,TV,TF,nim,varargin)
  % UNBAKE_NORMAL_MAP Given a uv-mapped (TV,TF) 3D surface (V,F) mesh and a
  % corresponding normal map image (nim), recover a high-resolution 3D geometry
  % that "explains" the normal map. Resolution (and runtime) is controlled by
  % the resolution of the normal map image (hint: start with a downsampled
  % image).
  %
  % [V,F,TV,TF,obj] = unbake_normal_map(V,F,TV,TF,nim,varargin)
  %
  % Inputs:
  %   V  #V by 3 list of 3D vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   TV  #TV by 3 list of 2D uv vertex positions
  %   TF  #F by 3 list of triangle indices into TV
  %   nim  #nim by #nim by 3 image of *object-space* normals (in range [0,1])
  % Outupts:
  %   V  #V by 3 list of 3D vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   TV  #TV by 3 list of 2D uv vertex positions
  %   TF  #F by 3 list of triangle indices into TV
  %
  % See also: normal_to_displacement_map, tangent_to_object_normal_map
  %

  function [sel,V,F,TV,TF] = upsample_helper(V,F,TV,TF)
    mean_area = mean(doublearea(TV,TF));
    pixel_area = 1./(size(nim,1)*size(nim,2));
    dblA = doublearea(TV,TF);
    sel = find(dblA>pixel_area);
    if isempty(sel)
      return;
    end
    % Selected edges on 3D mesh (V,F)
    [uE,~,EMAP] = unique(sort(reshape(F(:,[2 3 1 3 1 2]),[],2),2),'rows');
    uEM = sparse(EMAP,1,repmat(sparse(sel,1,1,size(F,1),1),3,1),size(uE,1),1)>0;
    MF = full(reshape(uEM(EMAP),[],3));
    % Selected edges on texture mesh (TV,TF)
    [uE,~,EMAP] = unique(sort(reshape(TF(:,[2 3 1 3 1 2]),[],2),2),'rows');
    uEM = sparse(EMAP,1,repmat(sparse(sel,1,1,size(TF,1),1),3,1),size(uE,1),1)>0;
    MTF = full(reshape(uEM(EMAP),[],3));
    M = MF | MTF;
    [TV,TF] = upsample(TV,TF,'OnlySelected',M);
    mprev = size(F,1);
    [V,F] = upsample(V,F,'OnlySelected',M);
    if ~quiet
      fprintf('#F: %g -> %g\n',mprev,size(F,1));
    end
  end

  max_iters = 1000;
  quiet = true;
  quiet = false;
  max_step_size_seen = -inf;
  method = 'lbfgs';
  progressive_upsampling = true;
  fix_symmetry = false;

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'FixSymmetry','Method','Quiet','Progressive'}, ...
    {'fix_symmetry','method','quiet','progressive_upsampling'});
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

  % Progressive Upsampling is better: similar amount of time for a far better
  % solution.
  if ~progressive_upsampling
    while true
      [sel,V,F,TV,TF] = upsample_helper(V,F,TV,TF);
      if isempty(sel)
        break;
      end
    end
  end

  cross2 = @(a,b,c) ...
    [a(:,2).*b(:,3)-a(:,3).*b(:,2), ...
     a(:,3).*b(:,1)-a(:,1).*b(:,3), ...
     a(:,1).*b(:,2)-a(:,2).*b(:,1)];
  nrmlev = @(p1,p2,p3) cross2(p2-p1,p3-p1);
  unrml = @(V,F) normalizerow(nrmlev(V(F(:,1),:),V(F(:,2),:),V(F(:,3),:)));


  while true
    % Compute normal averaged over all texture-space pixels landing in each
    % triangle: ~/uniform quadrature in texture space
    [X,Y] = meshgrid(linspace(0,1,size(nim,2)),linspace(1,0,size(nim,1)));
    I = in_element_aabb(TV(:,1:2),TF,[X(:) Y(:)]);
    Nnim = reshape(nim,[],3);
    TNsum = sparse( ...
      repmat(I(I>0),1,3),repmat(1:3,sum(I>0),1),Nnim(I>0,:),size(TF,1),3);
    TNcount =  sparse(I(I>0),1,1,size(TF,1),1);
    TN = normalizerow(full(TNsum./TNcount)*2-1);
    % Some triangles might not have received any samples
    N = unrml(V,F);
    TN(any(isnan(TN),2),:) = N(any(isnan(TN),2),:);
    A = diag(sparse(doublearea(V,F)));
    if fix_symmetry
      % flip if aggregate normal points in opposite direction from initial
      % normal
      D = A*sum(TN.*N,2);
      [C,CF] = connected_components(TF);
      CD = accumarray(CF',D)./accumarray(CF',diag(A));
      %tsurf(TF,[TV(:,1:2) C'],'CData',1*(CD(CF)<0.8));
      %T = normals(TV(:,1:2)*eye(2,3),TF);
      %tsurf(TF,[TV(:,1:2) C'],'CData',1*(T(:,3)>0));
      flipped = (CD(CF)<0.8);
      tsurf(TF(flipped,:),[TV(:,1:2) C'],'FaceVertexCData',N(flipped,:)*0.5+0.5);
      hold on;
      surf(X,Y,0*X,falpha(1,0),'FaceColor','texturemap','CData',nim);
      hold off;
      %TN(D<0,:) = -TN(D<0,:);
      error
    end
    ATN = A*TN;

    vec = @(X) X(:);
    Asqr = @(X) sum(sum(X.*(A*X)));
    obj_fun = @(V) Asqr(unrml(V,F)-TN);
    grad_fun = @(V) normal_gradient(V,F,vec(A*unrml(V,F))-vec(ATN));

    r3 = @(X) reshape(X,[],3);
    % This is awful...
    %options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','iter');
    %[V,obj1] = fminunc(@(V) deal(obj_fun(r3(V)),grad_fun(r3(V))),vec(V),options);

    switch method
    case 'lbfgs'
      % This thing actually seems to work
      [V,obj1] = fminlbfgs( ...
        @(V) fun_and_grad(V,@(V) obj_fun(r3(V)),@(V) grad_fun(r3(V))),vec(V), ...
        struct('Display','final','GradObj','on','TolFun',1e-6,'TolX',1e-5));
    case 'gd'
      %% My gradient descent is faster
      %[V,obj1] = fminlbfgs( ...
      %  @(V) fun_and_grad(V,@(V) obj_fun(r3(V)),@(V) grad_fun(r3(V))),vec(V), ...
      %  struct('HessUpdate','steepdesc','Display','final','GradObj','on','TolX',5e-5));
      [V,obj1] = fmingd( ...
        @(V) fun_and_grad(V,@(V) obj_fun(r3(V)),@(V) grad_fun(r3(V))),vec(V));
    end
    V = r3(V);

    %obj1 = obj_fun(V);
    %prev_step = 0;
    %step_size = 1;
    %prev_V = V;
    %for iter = 1:max_iters
    %  %method = 'gradient-descent';
    %  %method = 'stochastic-gradient-descent';
    %  %method = 'momentum';
    %  method = 'nesterov';
    %  %method = 'gauss-newton';
    %  %method = 'levenberg-marquardt';
    %  switch method
    %  case {'gradient-descent','gauss-newton','levenberg-marquardt','momentum'}
    %    % Recompute current normals and gradients
    %    dNdVT = normal_gradient(V,F);
    %    grad = dNdVT*(vec(A*unrml(V,F))-vec(ATN));
    %    switch method
    %    case {'gauss-newton','levenberg-marquardt'}
    %      max_step_size = 100;
    %      % This works really well for small meshes but takes forever/diverges for
    %      % big meshes
    %      K = (dNdVT*A*dNdVT');
    %      switch method 
    %      case 'levenberg-marquardt'
    %        K = (dNdVT*A*dNdVT');
    %        T = diag(diag(K));
    %        lambda = 1;
    %        K = K + lambda*T;
    %      end
    %      % This works really well for small meshes but takes forever/diverges for
    %      % big meshes
    %      step = reshape(K\(-grad),[],3);
    %    case 'gradient-descent'
    %      % Gradient descent
    %      step = -reshape(grad,[],3);
    %      max_step_size = 1;
    %    case 'momentum'
    %      % Fairly sensitive to mu
    %      mu = 0.5;
    %      step = mu*prev_step + (1-mu)*(-reshape(grad,[],3));
    %      % works better if I use the non line search step here...
    %      prev_step = step;
    %      max_step_size = 1;
    %    end
    %  case 'nesterov'
    %    % following AQP paper
    %    eta = 10;
    %    theta = (1-sqrt(1/eta))/(1+sqrt(1/eta));
    %    V = (1+theta)*V - theta*prev_V;
    %    dNdVT = normal_gradient(V,F);
    %    grad = dNdVT*(vec(A*unrml(V,F))-vec(ATN));
    %    step = -reshape(grad,[],3);
    %    max_step_size = 1;
    %  case 'stochastic-gradient-descent'
    %    % not impressed with this... maybe it'll be worthwhile if to keep things
    %    % in memory
    %    %
    %    % Especially doesn't seem to play well with the "stalling" criteria.
    %    batch_size = min(size(F,1),10000);
    %    P = randperm(size(F,1));
    %    P = P(1:batch_size);
    %    dNdVT = normal_gradient(V,F(P,:));
    %    grad = dNdVT*(vec(A(P,P)*unrml(V,F(P,:)))-vec(ATN(P,:)));
    %    step = -reshape(grad,[],3);
    %    max_step_size = 1;
    %  end
    %  prev_V = V;
    %  [step_size,V,obj,step] = line_search(obj_fun,@(V) V,V,step,max_step_size);
    %  max_step_size_seen = max(max_step_size_seen,step_size);
    %  V = reshape(V,[],3);
    %  if step_size == 0
    %    warning('step size == 0');
    %    break;
    %  end
    %  
    %  obj0 = obj1;
    %  obj1 = obj_fun(V);
    %  if ~quiet
    %    fprintf('  iter: % 5d, obj: %g\n',iter,obj1);
    %  end
    %  if obj0 - obj1 < 1e-5
    %    break
    %  end
    %end
    %if ~quiet
    %  fprintf('converged after %d iterations to %g\n',iter,obj1);
    %end


    if ~quiet
      tsurf(F,V,fsoft,'FaceVertexCData',repmat(orange,size(F,1),1),falpha(1,0.1));
      view(2);axis equal;camlight;
      %apply_ambient_occlusion([],'AddLights',false,'Factor',1);
      drawnow;
    end

    if progressive_upsampling
      [sel,V,F,TV,TF] = upsample_helper(V,F,TV,TF);
      if isempty(sel)
        break;
      end
    else
      % only one iteration needed: we upsampled everything at the beginning
      break;
    end
  end

  obj = obj_fun(V);
  if ~quiet
    fprintf('Final objective: %g\n',obj);
  end

  if ~quiet
    fprintf('max_step_size_seen: %g\n',max_step_size_seen);
  end

end
