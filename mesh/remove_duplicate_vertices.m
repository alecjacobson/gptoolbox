function [SV,SVI,SVJ] = remove_duplicate_vertices(V,epsilon,varargin)
  % REMOVE_DUPLICATE_VERTICES Remove duplicate vertices upto a uniqueness
  % tolerance (epsilon)
  %
  % SV = remove_duplicate_vertices(V,epsilon)
  % [SV,SVI,SVJ] = ...
  %   remove_duplicate_vertices(V,epsilon,'ParameterName',ParameterValue,...)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   epsilon  uniqueness tolerance (significant digit)
  %   Optional:
  %     'WhiteList' Only merge vertices from the following selection (not
  %     working correctly, yet)
  % Outputs:
  %   SV  #SV by dim new list of vertex positions
  %   SVI #SV by 1 list of indices so SV = V(SVI,:) 
  %   SVJ #V by 1 list of indices so V = SV(SVJ,:)
  %
  % Example:
  %   % Mesh in (V,F)
  %   [SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-7);
  %   % remap faces
  %   SF = SVJ(F);
  %

  % default values
  whitelist = [];
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'WhiteList'}, ...
    {'whitelist'});
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

  if isempty(whitelist)
    assert(nargin==1 || epsilon >= 0);
    if nargin==1 || epsilon == 0
      [SV,SVI,SVJ] = unique(V,'rows','stable');
    else
      [~,SVI,SVJ] = unique(round(V/(10*epsilon)),'rows','stable');
      SV = V(SVI,:);
    end
  else
    error('not implemented correctly')
    VW = V(whitelist,:);
    J = 1:size(V,1);
    JW = J(whitelist);
    [SVW,SVIW,SVJW] = remove_duplicate_vertices(VW,epsilon);
    SJ = 1:size(V,1);
    SJ(whitelist) = JW(SVJ);
  end
end
