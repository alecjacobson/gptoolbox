function [NN,EE,A] = spherical_subdivision(N,E,varargin)
  % SPHERICAL_SUBDIVISION  Given a "spherical edge network", subdivide the
  % curves while keeping all points on the (unit) sphere.
  %
  % [NN,EE] = spherical_subdivision(N,E);
  % [NN,EE] = spherical_subdivision(N,E,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   N  #N by 3 Vertices on the unit sphere
  %   E  #E by 2 list of edge indices into N
  %   Optional:
  %     'MinLength' followed by minimum edge length to subdivide {0}.
  %     'MaxIter' followed by maximum number of iterations {1}.
  % Outputs:
  %   NN  #NN by 3 Vertices on the unit sphere
  %   EE  #E by 2 list of edge indices into N
  %   A  #NN by #N subdivision matrix
  %
  
  % default values
  min_length = 0;
  max_iter = 1;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'MinLength','MaxIter',}, ...
    {'min_length','max_iter'});
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

  NN = N;
  EE = E;
  A = speye(size(N,1),size(N,1));
  iter = 1;
  while true
    l = sqrt(sum((NN(EE(:,1),:) - NN(EE(:,2),:)).^2,2));
    long = l>min_length;
    newNI = size(NN,1)+(1:sum(long))';
    NN = [NN;slerp(NN(EE(long,1),:),NN(EE(long,2),:),0.5)];
    A(newNI,:) = 0.5*(A(EE(long,1),:) + A(EE(long,2),:));
    newEE = [EE(long,1) newNI;newNI EE(long,2)];
    EE = [EE(~long,:);newEE];
    if iter >= max_iter
      break;
    end
    iter = iter + 1;
  end

end
