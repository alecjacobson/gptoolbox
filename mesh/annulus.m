function [V,F,bi,bo] = annulus(s,r,varargin)
  % ANNULUS Construct a triangle mesh of a unit annulus.
  % 
  % [V,F] = annulus(s,r)
  % [V,F] = annulus(s,r,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   s  number of samples on the inner ring
  %   r  radius of the inner ring
  %   Optional:
  %     'Flags'  followed by flags to pass to Triangle for meshing. 
  %        {'-q30 -aX'} where X is squared boundary edge length
  %     'R'  followed by outer ring radius {1}
  % Outputs:
  %   V  #V by 2 list of mesh vertex positions
  %   F  #F by 3 list of triangle mesh indices
  %   bi  list of vertices on inner boundary
  %   bo  list of vertices on outer boundary
  %   

  flags = nan;
  R = 1;
  params_to_variables = containers.Map( ...
    {'Flags','R'},{'flags','R'});
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
  if isnan(flags)
    flags = sprintf('-q30 -a%0.17f',(2*pi*r/s)^2);
  end

  theta = linspace(0,2*pi,s+1)';theta = theta(1:end-1);
  Vr = r*[cos(theta) sin(theta)];
  Er = fliplr([1:size(Vr,1);2:size(Vr,1) 1]');
  theta = linspace(0,2*pi,ceil(R/r)*s+1)';theta = theta(1:end-1);
  VR = R*[cos(theta) sin(theta)];
  ER = [1:size(VR,1);2:size(VR,1) 1]';

  [V,F] = triangle([Vr;VR],[Er;size(Vr,1)+ER],[0 0],'Flags',flags);
  b = unique(outline(F));
  bi = intersect(find(normrow(V)< 0.5*(R+r)),b);
  bo = intersect(find(normrow(V)> 0.5*(R+r)),b);
end
