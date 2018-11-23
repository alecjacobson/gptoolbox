function [SE,II] = sharp_edges(V,F,varargin)
  % SHARP_EDGES Given a mesh, compute sharp edges.
  %
  % [SE] = sharp_edges(V,F)
  % [SE,II] = sharp_edges(V,F,'ParameterName',ParameterValue, ...)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions 
  %   F  #F by 3 list of mesh triangle indices into V
  %   Optional:
  %     'Angle'  followed by dihedral angle considered sharp. {pi*0.11}
  % Outputs:
  %   SE  #SE by 2 list of edge indices into V 
  %   II  #II by 2 list of indices into F so that: 
  %     SE=unique(sort(F(II),2),'rows');
  %

  angle = pi*0.11;
  % default values
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Angle'},{'angle'});
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

  % sharp edges
  [A,C] = adjacency_dihedral_angle_matrix(V,F);
  %% This is much much slower
  %A(1&A) = abs(A(1&A)-pi)>pi*0.11;
  [AI,AJ,AV] = find(A);
  keep = abs(AV-pi)>(angle) & ~isnan(AV);
  A = sparse(AI(keep),AJ(keep),1,size(A,1),size(A,2));
  [CI,~,CV] = find(C.*A);
  II = [CI+mod(CV,3)*size(F,1) CI+mod(CV+1,3)*size(F,1)];
  E = F(II);
  SE = unique(sort(E,2),'rows');

end
