function C = catmull_rom_eval(P,T,t,varargin)
  % C = catmull_rom_eval(P,T,t,varargin)
  %
  % Inputs:
  %   P  #P by dim list of control points
  %   T  #T list of control point parametric values (times)
  %   t  #t list of parametric evaluation values (times)
  %   Optional:
  %     'Period'  followed by finite period for closed curve {inf → open}
  %     'Tao'  followed by tension parameter {0.5}
  % Outputs:
  %   C  #t by dim list of evaluation points
  %

  tao = 0.5;
  period = inf;
  params_to_variables = containers.Map( ...
    {'Period','Tao'}, ...
    {'period','tao'});
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


  n = size(P,1);
  % https://github.com/alecjacobson/computer-graphics-kinematics-solution/blob/master/catmull_rom_interpolation.cpp
  T = reshape(T,[],1);
  assert(n == size(T,1));
  assert(all(t>=min(T)));
  if period==inf
    assert(all(t<=max(T)));
  else
    assert(min(t) >= 0);
    t = mod(t,period);
  end
  % find keyframe just before s
  [~,K1] = histc(t,[T;period+T(1)]);
  if period==inf
    assert(all(K1~=0));
    K1(t==max(T)) = n-1;
  end
  K2 = mod((K1+1)-1,n)+1;

  % parametric distance between k₁ and k₂
  s = (t-T(K1))./mod(T(K2) - T(K1),period);
  % values at k₁ and k₂  
  P1 = P(K1,:);
  P2 = P(K2,:);
  if period == inf
    D1 = tao*0.5*(P2-P1);
    D1(K1>1,:) = tao*(P2(K1>1,:) - P(K1(K1>1)-1,:));
    D2 = tao*0.5*(P2-P1);
    D2(K2<n,:) = tao*(P(K2(K2<n)+1,:) - P1(K2<n,:));
  else
    K0 = mod(K1-1 - 1,n)+1;
    K3 = mod(K2+1 - 1,n)+1;
    D1 = tao*(P2 - P(K0,:));
    D2 = tao*(P(K3,:) - P1);
  end
  M = [0 1 0 3;0 1 0 2;0 1 1 1;1 1 0 0]';
  C = sum(permute([s.^(3:-1:0)]/M,[1 3 2]).*cat(3,P1,P2,D1,D2),3);

end
