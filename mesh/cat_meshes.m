function [V,F] = cat_meshes(V,F,varargin)
  % CAT_MESHES concatenate many meshes
  %
  % [V,F] = cat_meshes(V1,F1,V2,F2, ...
  % 
  % Inputs:
  %   V1  #V by dim list of vertex positions of mesh 1
  %   F1  #F by simplex-size list of simplex indices of mesh 1
  %   V2  #V by dim list of vertex positions of mesh 2
  %   F2  #F by simplex-size list of simplex indices of mesh 2
  %   ...
  % Outputs:
  %   V  #V1+#V2+... by dim list of vertex positions
  %   F  #F1+#F2+... by simplex-size list of simplex indices
  %
  assert(mod(nargin,2)==0);
  for m = 1:nargin/2
    Vm = varargin{2*m-1};
    Fm = varargin{2*m};
    F = [F; Fm+size(V,1)];
    V = [V; Vm];
  end
end

