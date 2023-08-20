function [V,F] = buckliball(varargin)
  % [V,F] = buckliball();
  % [V,F] = buckliball('ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   Optional:
  %      'h'  followed by the desired edge length {0.1}
  %      'theta'  followed by the desired hole angle (0,π/10) {pi*0.9}
  %      'Vis' followed by whether or not to visualize {false}
  %      'MaxIter' followed by the maximum number of iterations {100}
  %      'Tol' followed by the tolerance for lloyd relaxation convergence {0.1*h}
  % Outputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle indices into V
  %


  vis = false;
  h = 0.1;
  theta = pi * 0.09;
  tol = 0.1*h;
  max_iter = 100;
  params_to_variables = containers.Map( ...
    {'h','theta','Vis','MaxIter','Tol'}, ...
    {'h','theta','vis','max_iter','tol'});
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

  if theta>pi*0.1
    warning('Theta is larger than pi*0.1. Holes will overlap');
  end

  [V,F,I,C] = rhombicosidodecahedron();
  Q = find(accumarray(I,1)==2);
  N = C(Q,:);

  CV = [];
  CE = [];
  for i = 1:size(N,1)
    [iV,iE] = circle_of_latitude(N(i,:),theta,h);
    CE = [CE size(CV,2)+iE'];
    CV = [CV iV'];
  end
  CV = CV';
  CE = CE';


  if vis
    tsurf(CE,CV);
    axis equal;
    drawnow;
    error
  end

  rn = ceil( 1.0 * ( 5*pi*1/h^2-size(CV,1) ) );

  % Better initial sampling
  ssn = ceil(-(log(3)-log(3/10*rn-3/5))/log(3));
  [SV,SF] = subdivided_sphere(ssn, 'SubdivisionMethod','sqrt3');
  V = SV;
  %V = randsphere( rn );


  V = [CV;V];
  V0 = V;
  b = (1:size(CV,1))';

  for iter = 1:max_iter
    % much faster than convhull. Why?
    F = convhulln(V);

    if vis
      tsurf(F,V);
      hold on;
      sct(V(b,:),'or','LineWidth',2);
      hold off;
      axis equal;
      drawnow
    end

    BC = barycenter(V,F);
    A = doublearea(V,F);

    B = sparse(F,repmat(1:size(F,1),3,1)',[A A A],size(V,1),size(F,1));
    V0 = V;
    V = normalizerow( (B*BC)./(B*ones(size(F,1),1)) );
    V(b,:) = V0(b,:);
    small_change = norm(V - V0,inf) < tol;

    % Flip edges crossing hole boundaries
    while true
      W = sign(winding_number([CV;0 0 0],[repmat(size(CV,1)+1,size(CE,1),1) CE],V,'Fast',true));
      E = edges(F);
      X = W(E(:,2))~=W(E(:,1));
      X(any(ismember(E,b),2)) = false;
      if ~any(X)
        break;
      end
      F = flip_edges(F,E(X,:));
    end

    % Remove triangles inside each hole
    BC = barycenter(V,F);
    W = winding_number([CV;0 0 0],[repmat(size(CV,1)+1,size(CE,1),1) CE],BC,'Fast',true);
    F = F(W>0,:);
    [V,I,~,F] = remove_unreferenced(V,F);
    b = I(b);

    if small_change
      break;
    end

    if iter == max_iter
      warning('Max iterations (%d) reached',max_iter);
    end

    %tsurf(F,V,'CData',doublearea(V,F));
    %axis equal;
    %drawnow;

  end


  if vis
    tsurf(F,V);
    hold on;
    sct(V(b,:),'or','LineWidth',2);
    hold off;
    axis equal;
  end

  function [V,E] = circle_of_latitude(N,theta,h)
    r = sin(theta);
    % circumferential length
    n = ceil(2*pi*r/h);
    assert(n,'h too small');
    phi = linspace(0,2*pi,n+1)';
    phi = phi(1:end-1);
    V = [r*[cos(phi) sin(phi)] repmat(cos(theta),n,1)];
    E = [1:n;2:n 1]';
    % rotate so that [0 0 1] points to N
    R = R_from_N(N);
    V = V*R;
  end

  function R = R_from_N(N)
    if N(3) == 1
      R = eye(3);
    elseif N(3) == -1
      R = [1 0 0;0 -1 0;0 0 -1];
    else
      [~,~,R] = axisanglebetween(N,[0 0 1]);
    end
  end

end
