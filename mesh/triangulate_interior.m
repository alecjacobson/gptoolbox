function [TV,TF] = triangulate_interior(V,E,varargin)
  % TRIANGULATE_INTERIOR  Triangulate the interior of a polygon or set of
  % polygons.
  %
  % [TV,TF] = triangulate_interior(V,E,varargin)
  %
  % Inputs:
  %   V  #V by 2 list of vertex positions
  %   E  #E by 2 list of edge indices into V
  %   Optional:
  %     'Flags'  followed by flags to pass to refine_triangulation
  %     'Filter'  followed by a function handle to filter out faces of initial
  %       pass based on winding number with respect to (V,E) computed at
  %       barycenters. {@(W) abs(W)>0.5}
  % Outputs:
  %   TV #TV by 2 list of vertex positions
  %   TF #TF by 3 list of triangle indices into TV
  %
  % See also: triangulate, refine_triangulation
  %


  flags = '';
  filter = @(W) abs(W)>0.5;
  BI = (1:size(E,1))';
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'Flags','Filter','BoundaryEdges'}, ...
    {'flags','filter','BI'});
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

  [V,~,~,E] = remove_duplicate_vertices(V,0,'F',E);
  E = E(E(:,1)~=E(:,2),:);

  % Sometimes meshing the convex hull can add super slivers which are
  % numerically annoying to remove with the winding number. Avoid this if we
  % can (i.e., the input already defines a closed outer hull.
  [Vc,Fc] = triangulate(V,E,'Flags','');
  if isempty(Fc) || numel(unique(Fc)) < size(V,1)
    [Vc,Fc] = triangulate(V,E,'Flags','c');
  end
  BC = barycenter(Vc,Fc);
  W = winding_number(V,E(BI,:),BC);
  Fc = Fc(filter(W),:);

  [Vc,~,~,Fc] = remove_unreferenced(Vc,Fc,true);

  %tsurf(Fc,Vc,'FaceColor',blue);
  %pause
  [TV,TF] = refine_triangulation(Vc,E,Fc,'Flags',flags);

end
