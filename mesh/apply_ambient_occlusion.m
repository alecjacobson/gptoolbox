function [AO,l] = apply_ambient_occlusion(t,varargin)
  % APPLY_AMBIENT_OCCLUSION  Apply ambient occlusion to a given plot of a mesh.
  % Derive colors from the current plot and provide new modulated colors
  % according to the ambient occlusion computed for this mesh.
  %
  % AO = apply_ambient_occlusion(t)
  % [AO,l] = apply_ambient_occlusion(t,'ParameterName',parameter_value, ...)
  %
  % Inputs:
  %   t  handle to a `trisurf` or `patch` object
  %   Optional:
  %     'AO' followed by previously computed ambient occlusion values {[]}
  %     'AddLights' followed by whether to add lights and soft lighting {true}
  %     'Samples' followed by the number of samples to use when computing
  %       ambient occlusion {1000}
  % Outputs:
  %   AO  #V by 1 list of ambient occlusion values
  %   l  #l list of light handles
  AO = [];
  % default values
  add_lights = true;
  samples = 1000;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'AO','AddLights','Samples'}, ...
    {'AO','add_lights','samples'});
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
  V = t.Vertices;
  Poly = t.Faces;
  % triangulate high order facets
  F = [];
  for c = 3:size(Poly,2)
    F = [F;Poly(:,[1 c-1 c])];
  end
  C = t.FaceVertexCData;
  if size(C,2) == 1
    C = squeeze(ind2rgb(floor(matrixnormalize(t.FaceVertexCData)*size(colormap,1))+1,colormap));
  end
  if size(C,1) == size(V,1)
    t.FaceColor = 'interp';
    O = V;
    N = per_vertex_normals(V,F);
    % Matlab uses backwards normals
    t.VertexNormals =  -N;
    lighting phong;
  elseif size(C,1) == size(F,1)
    t.FaceColor = 'flat';
    % Could be fancy and use high-order quadrature.
    O = barycenter(V,F);
    N = normalizerow(normals(V,F));
    lighting flat;
  else
    error('Unknown number of colors');
  end
  if isempty(AO)
    AO = ambient_occlusion(V,F,O,N,samples);
  end
  t.FaceVertexCData = bsxfun(@times,C,1-AO);
  if add_lights
    t.FaceLighting = 'phong';
    l = {};
    l{1} = light('Position',[1 1 100],'Style','infinite');
    l{2} = light('Position',[1 100 1],'Style','infinite');
    t.SpecularStrength = 0.1;
    t.DiffuseStrength = 0.1;
    t.AmbientStrength = 0.8;
  end
end
