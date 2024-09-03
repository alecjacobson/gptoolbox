function [err,err1_2,err2_1] = spline_spline_chamfer(P1,C1,P2,C2,varargin)
  % SPLINE_SPLINE_CHAMFER Compute chamfer distance between two splines.
  % Specifically this is an approximation of the integrated sum of squared
  % distances.
  % 
  % err1_2 = ∫_S₁ min_{x₂ ∈ S₂} ||x₁ - x₂||² dx₁ / ∫_S₁ dx₁
  %
  % Note that the integral is in the geometric space of the spline, not the
  % parameter space.
  %
  % [err,err1_2,err2_1] = spline_spline_chamfer(P1,C1,P2,C2,varargin)
  %
  % Inputs:
  %  P1  #P1 by dim list of control points
  %  C1  #C1 by 4 list of cubic Bezier control points indices into P1
  %  P2  #P2 by dim list of control points
  %  C2  #C2 by 4 list of cubic Bezier control points indices into P2
  %  Optional:
  %    'Method' followed by one of:
  %      'gauss-legendre'  use gauss-legendre quadrature to sample the spline
  %      'uniform'  use uniform spline sampling to sample
  %      'monte-carlo'  use monte-carlo sampling on a piece-wise linear
  %        approximation of the spline
  %    'Samples' followed by number of samples
  %    'Tol' followed by flatness tolerance
  % Outputs:
  %   err  total chamfer distance
  %   err1_2  chamfer distance from spline 1 to spline 2
  %   err2_1  chamfer distance from spline 2 to spline 1
  %

  function [A,a] = spline_arc_lengths(P,C)
    a = 0;
    nq = 10;
    A = zeros(size(C,1),1);
    for c = 1:size(C,1)
      A(c) = cubic_arc_length(P(C(c,:),:),nq,0,1);
      a = a + A(c);
    end
  end

  function [U1,W1,a1] = spline_gauss_legendre(P1,C1,num_samples)
    [A1,a1] = spline_arc_lengths(P1,C1);
    S1 = round(num_samples.*A1./a1);
    U1 = [];
    W1 = [];
    for c = 1:size(C1,1)
      [Tc,Wc] = gauss_legendre_quadrature(S1(c),0,1);
      Uc = cubic_eval(P1(C1(c,:),:),Tc);
      U1 = [U1; Uc];
      W1 = [W1; A1(c)*Wc];
    end
  end

  function [U1,W1,a1] = spline_uniformly_sample(P1,C1,num_samples)
    [A1,a1] = spline_arc_lengths(P1,C1);
    S1 = round(num_samples.*A1./a1);
    U1 = [];
    W1 = [];
    for c = 1:size(C1,1)
      %[~,Fc,Uc] = cubic_uniformly_sample(P1(C1(c,:),:),2*S1(c)+1);
      %Fc = Fc(1:2:end) + Fc(2:2:end);
      %Uc = Uc(2:2:end,:);
      %U1 = [U1; Uc];
      %W1 = [W1; Fc];

      % This seems to work better/same for small number of samples
      [~,Fc,Uc] = cubic_uniformly_sample(P1(C1(c,:),:),S1(c)+1);
      Fc = Fc * 0.5;
      Wc = [Fc;0] + [0;Fc];
      U1 = [U1; Uc];
      W1 = [W1; Wc];
    end
    % This isn't even worth it.
    %[U1,~,I] = unique(U1,'rows');
    %W1 = accumarray(I,W1,[size(U1,1),1]);
  end

  %num_samples = 10000;
  %method = 'monte-carlo';
  num_samples = 100;
  method = 'uniform';
  flat_tol = 1e-4;
  symmetric = true;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Symmetric','Method','Samples','Tol'}, ...
    {'symmetric','method','num_samples','tol'});
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

  %% Dense Sample
  if symmetric || strcmp(method,'monte-carlo')
    [V1,E1] = spline_to_poly(P1,C1,flat_tol);
  end
  [V2,E2] = spline_to_poly(P2,C2,flat_tol);
  switch method
  case 'gauss-legendre'
    % Not impressed with this here. Seems to be running into numerical issues
    % before it converges?
    [U1,W1,a1] = spline_gauss_legendre(P1,C1,num_samples);
    if symmetric
      [U2,W2,a2] = spline_gauss_legendre(P2,C2,num_samples);
    end
  case 'uniform'
    % This works well for relatively small number of samples
    [U1,W1,a1] = spline_uniformly_sample(P1,C1,num_samples);
    if symmetric
      [U2,W2,a2] = spline_uniformly_sample(P2,C2,num_samples);
    end
  case 'monte-carlo'
    % Fast but requires a lot of samples
    % num_samples = 100000;
    %num_samples1 = round(max(size(E1,1),num_samples));
    %num_samples2 = round(max(size(E2,1),num_samples));
    
    num_samples1 = num_samples;
    num_samples2 = num_samples;
    
    U1 = rand_samples(V1,E1,num_samples1);
    a1 =  sum(edge_lengths(V1,E1));
    W1 = repmat(a1/num_samples1,num_samples1,1);

    if symmetric
      U2 = rand_samples(V2,E2,num_samples2);
      a2 =  sum(edge_lengths(V2,E2));
      W2 = repmat(a2/num_samples2,num_samples2,1);
    end
    
  end

  %% Compute Chamfer Distance
  squared_dist1_2 = point_mesh_squared_distance(U1,V2,E2,'Method','libigl');
  %% Why is this not reliable?
  %squared_dist2_1 = point_spline_squared_distance(U2,P1,C1,flat_tol);
  %squared_dist1_2 = point_spline_squared_distance(U1,P2,C2,flat_tol);

  %err = 0.5*sum(squared_dist2_1) / sum(edge_lengths(V2,E2)) + ...
        %0.5*sum(squared_dist1_2) / sum(edge_lengths(V1,E1));
  err1_2 = 0.5*sum(squared_dist1_2.*W1) / a1;
  if symmetric
    squared_dist2_1 = point_mesh_squared_distance(U2,V1,E1,'Method','libigl');
    err2_1 = 0.5*sum(squared_dist2_1.*W2) / a2;
    err = err1_2 + err2_1;
  else
    err = err1_2;
    err2_1 = nan;
  end
end
