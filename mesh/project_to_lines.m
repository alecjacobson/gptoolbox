function [T,sqrD] = project_to_lines(P,S,D,varargin)
  % PROJECT_TO_LINES  project points onto vectors, that is find the paramter
  % t for a point p such that proj_p = (y-x).*t, additionally compute the
  % sqaured distance from p to the line of the vector, such that 
  % |p - proj_p|² = sqr_d
  %
  % [T,sqrD] = project_to_lines(P,S,D)
  % [T,sqrD] = project_to_lines(P,S,D,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   P  #P by dim list of points to be projected
  %   S  #vectors by dim list of start positions of each vector
  %   D  #vectors by dim list of destination positions of each vector
  %   Optional:
  %     'Segments' followed by whether to project onto line segments (between S
  %       and D) {false}
  % Outputs:
  %   T  #P by #vectors list of parameters for each pair of points in p and
  %     vectors in (S,D)
  %   sqrD  #P by #vectors list of squared distances for each paint of points
  %     in p and vectors in (S,D)
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %

  % default values
  segments = false;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Segments'}, ...
    {'segments'});
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

  % number of dimensions
  dim = size(P,2);
  assert(dim == size(S,2));
  assert(dim == size(D,2));
  % number of vectors
  nv = size(S,1);
  assert(nv == size(D,1),'#S should equal #D');
  % number of points
  np = size(P,1);
  if(nv == 1 && np == 1)
    % optimize case where we're just handling one point projected onto one line
    DmS = D-S;
    v_sqrlen = sum((DmS).^2,2);
    T = - sum((S-P).*(DmS)) ./ v_sqrlen;
    if segments
      T = min(max(T,0),1);
    end
    if(nargout > 1)
      sqrD = sum((P - ( (1-T) .* S + T .* D)).^2);
    end
  else
  
    % http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
    % points as seen by each edge P(i,j,:) gives coordinates of
    % mesh vertex i seen by jth vector
    PV = permute(repmat(P,[1 1 nv]),[1 3 2]);
    % Closest vertex positions of edge start and end points on domain mesh as 
    % seen by each domain vertex E1P(i,j,:) gives coordinates of
    % start point coordinates of jth edge as seen by mesh vertex  i
    SP = permute(repmat(S,[1 1 np]),[3 1 2]);
    DP = permute(repmat(D,[1 1 np]),[3 1 2]);
    % sqred length of each edge, v_sqrlen(j) gives squared length of each edge
    v_sqrlen = sum((S-D).^2,2);
    % Parameter t for each point on each edge's line equation, t(i,j) tells where
    % to find projection of mesh vertex j onto line of edge i, parameterized
    % linearly from start point (t=0) to end point (t=1)
    T = - dot(SP-PV,DP-SP,3)./ repmat(v_sqrlen',[np 1]);

    if segments
      T = min(max(T,0),1);
    end
    if(nargout > 1)
      % t seen by each coordinate
      Tdim = repmat(T,[1 1 dim]);
      % sqrD(i,j) tells squared distance from point i to line of vector j
      sqrD = sum((PV - ( (1-Tdim) .* SP + Tdim .* DP)).^2,3);
    end
  end
end
