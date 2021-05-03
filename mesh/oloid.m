function [V,F] = oloid(n,varargin)
  % OLOID a minimal mesh of an oloid (constructed via the convex hull of two
  % unit-radius circles
  %
  % Inputs:
  %   n  number of vertices sampled on each circle
  %   Optional:
  %     'Spacing' followed by distance between circle centers {0.5}
  %        0 → sphericon
  %        1 → oloid
  %        sqrt(2) → "two circle roller" with constant centroid height during
  %          rolling
  % Outputs:
  %   V  #V by 3 list of vertex positions (#V ≈ 1⅓ n)
  %   F  #F by 3 list of triangle indices into V (#F ≈ 2⅔ n)
  %
  h = 1;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Spacing'}, ...
    {'h'});
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

  th = linspace(0,2*pi,n+1)';
  th = th(1:end-1);

  VL = [cos(th)-h*0.5 sin(th) 0*th];
  VR = [cos(th)+h*0.5 0*th sin(th)];
  VL = VL(VL(:,1)<=0,:);
  VR = VR(VR(:,1)>=0,:);
  V = [VL;VR];
  % F = convhull(V);:WAY SLOWER!!
  F = convhulln(V);
  [V,~,~,F] = remove_unreferenced(V,F);

end
