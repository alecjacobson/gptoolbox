function [N,cN,cF] = per_corner_normals(V,F,varargin)
  % N = per_corner_normals(V,F,varargin)
  %
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions 
  %   F  #F by 3 list of mesh triangle indices into V
  %   Optional:
  %     'Angle'  followed by dihedral angle considered sharp. {pi*0.11}
  % Outputs:
  %   N  3*#F by 3 list of per-corner normals
  %   cN  #cN by 3 list of per-corner normals before expanding to corners
  %   cF  #F by 3 list of indices into cN
  % 
  % Example:
  %   [V,F] = cylinder_mesh(1,20,'Caps',true);
  %   N = per_corner_normals(V,F);
  %   tsurf(F,V,'FaceColor','w',falpha(1,0.5));
  %   hold on;
  %   qvr(V(F,:),N);
  %   hold off;
  %   axis equal;
  %
  % Example:
  %   [V,F] = cylinder_mesh(1,20,'Caps',true);
  %   N = per_corner_normals(V,F);
  %   % zip into gltf style mesh (shared F, per-vertex attributes)
  %   [uVN,~,I] = unique([V(F,:) N],'rows');
  %   uV = uVN(:,1:3);
  %   uN = uVN(:,3+(1:3));
  %   uF = reshape(1:size(F,1)*3,size(F));
  %   uF = I(uF);
  %
  % Example:
  %   [V,F] = cylinder_mesh(1,20,'Caps',true);
  %   % zip into an obj style mesh (indexed V,N with separate faces)
  %   [~,cN,cF] = per_corner_normals(V,F);
  %   writeOBJ('per_corner_normals-test.obj',V,F,[],[],cN,cF);
  %   
  % 

  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Angle'}, ...
    {'angle'});
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

  E = sharp_edges(V,F);
  [cF,I] = cut_edges(F,E);
  cV = V(I,:);
  cN = per_vertex_normals(cV,cF);
  N = cN(cF,:);
end
