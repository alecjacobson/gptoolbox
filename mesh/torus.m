function [V,F,Q] = torus(n,m,r,varargin)
  % TORUS Construct a triangle mesh of a unit torus.
  % 
  % [V,F] = torus(n,m,r)
  % [V,F] = torus(n,m,r,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   n  number of vertices around inner ring
  %   m  number of vertices around outer ring
  %   r  radius of the inner ring
  %   Optional:
  %     'R'  followed by outer ring radius {1}
  % Outputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of triangle mesh indices
  %
  % Example:
  %   % Roughly even shaped triangles
  %   n = 40;
  %   r = 0.4;
  %   [V,F] = torus(n,round(r*n),r);
  R = 1;
  params_to_variables = containers.Map( ...
    {'R'},{'R'});
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

  [V,F] = create_regular_grid(n,m,true,true);

  V = V*2*pi;
  th = V(:,2);
  phi = V(:,1);
  V = [cos(phi).*(R+r*cos(th)) sin(phi).*(R+r*cos(th)) r*sin(th)];
  Q = [F(1:2:end-1,[1 2]) F(2:2:end,[2 3])];


end
