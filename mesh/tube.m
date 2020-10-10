function [V,F] = tube(s,r,varargin)
  % TUBE Construct a triangle mesh of a unit tube.
  % 
  % [V,F] = tube(s,r)
  % [V,F] = tube(s,r,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   s  number of samples on the inner ring
  %   r  radius of the inner ring
  %   Optional:
  %     'Flags'  followed by flags to pass to Triangle for meshing. 
  %        {''} where X is squared boundary edge length
  %     'R'  followed by outer ring radius {1}
  %     'Levels'  followed by number of "stacks" or "levels" {1}
  % Outputs:
  %   V  #V by 2 list of mesh vertex positions
  %   F  #F by 3 list of triangle mesh indices
  %   

  flags = nan;
  R = 1;
  levels = 1;
  params_to_variables = containers.Map( ...
    {'Flags','R','Levels'},{'flags','R','levels'});
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
    %flags = sprintf('-q30 -a%0.17f',(2*pi*r/s)^2);
    flags = '';
  end
  [AV,AF] = annulus(s,r,'Flags',flags,'R',R);
  [V,F] = extrude(AV,AF,'Levels',levels);

end
