function writeDAE(filename,varargin)
  % WRITEDAE  Write a mesh to a Collada .dae scene file.
  %
  % writeDAE(filename,V,F)
  % writeDAE(filename,V1,F1,V2,F2, ...)
  %
  % Inputs:
  %   filename  path to .dae file
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %

  function s = id(p,n,i)
    s = sprintf('%sID%d',p,n+(i-1)/2*10);
  end

  docNode = com.mathworks.xml.XMLUtils.createDocument('COLLADA');
  docRootNode = docNode.getDocumentElement;
  docRootNode.setAttribute('xmlns','http://www.collada.org/2005/11/COLLADASchema');
  docRootNode.setAttribute('version','1.4.1');
  asset = docNode.createElement('asset');
  unit = docNode.createElement('unit');
  unit.setAttribute('meter','0.0254000');
  unit.setAttribute('name','inch');
  asset.appendChild(unit);
  up_axis = docNode.createElement('up_axis');
  up_axis.appendChild(docNode.createTextNode('Y_UP'));
  asset.appendChild(up_axis);
  docRootNode.appendChild(asset);

  visual_scenes = docNode.createElement('library_visual_scenes');
  visual_scene = docNode.createElement('visual_scene');
  visual_scene.setAttribute('id','ID2');
  sketchup = docNode.createElement('node');
  sketchup.setAttribute('name','SketchUp');

  library_nodes = docNode.createElement('library_nodes');

  library_geometries = docNode.createElement('library_geometries');

  for i = 1:2:numel(varargin)
    V = varargin{i};
    F = varargin{i+1};

    node = docNode.createElement('node');
    node.setAttribute('id',id('',3,i));
    node.setAttribute('name',sprintf('instance_%d',i-1));

    matrix = docNode.createElement('matrix');
    matrix.appendChild(docNode.createTextNode(sprintf('%d ',eye(4))));
    node.appendChild(matrix);
    instance_node = docNode.createElement('instance_node');
    instance_node.setAttribute('url',id('#',4,i));
    node.appendChild(instance_node);
    sketchup.appendChild(node);

    node = docNode.createElement('node');
    node.setAttribute('id',id('',4,i));
    node.setAttribute('name',sprintf('skp%d',i-1));
    instance_geometry = docNode.createElement('instance_geometry');
    instance_geometry.setAttribute('url',id('#',5,i));
    bind_material = docNode.createElement('bind_material');
    bind_material.appendChild(docNode.createElement('technique_common'));
    instance_geometry.appendChild(bind_material);
    node.appendChild(instance_geometry);
    library_nodes.appendChild(node);

    geometry = docNode.createElement('geometry');
    geometry.setAttribute('id',id('',5,i));
    mesh = docNode.createElement('mesh');
    source = docNode.createElement('source');
    source.setAttribute('id',id('',6,i));
    float_array = docNode.createElement('float_array');
    float_array.setAttribute('id',id('',7,i));
    float_array.setAttribute('count',num2str(numel(V)));
    float_array.appendChild(docNode.createTextNode(sprintf('%g ',V')));
    source.appendChild(float_array);
    technique_common = docNode.createElement('technique_common');
    accessor = docNode.createElement('accessor');
    accessor.setAttribute('count',num2str(size(V,1)));
    accessor.setAttribute('source',id('#',7,i));
    accessor.setAttribute('stride','3');
    for name = {'X','Y','Z'}
      param = docNode.createElement('param');
      param.setAttribute('name',name);
      param.setAttribute('type','float');
      accessor.appendChild(param);
    end
    technique_common.appendChild(accessor);
    source.appendChild(technique_common);
    mesh.appendChild(source);
    vertices = docNode.createElement('vertices');
    vertices.setAttribute('id',id('',8,i));
    input = docNode.createElement('input');
    input.setAttribute('semantic','POSITION');
    input.setAttribute('source',id('#',6,i));
    vertices.appendChild(input);
    mesh.appendChild(vertices);
    triangles = docNode.createElement('triangles');
    triangles.setAttribute('count',num2str(size(F,1)));
    input = docNode.createElement('input');
    input.setAttribute('offset','0');
    input.setAttribute('semantic','VERTEX');
    input.setAttribute('source',id('#',8,i));
    triangles.appendChild(input);
    p = docNode.createElement('p');
    p.appendChild(docNode.createTextNode(sprintf('%d ',(F-1)')));
    triangles.appendChild(p);
    mesh.appendChild(triangles);
    geometry.appendChild(mesh);
    library_geometries.appendChild(geometry);
  end

  visual_scene.appendChild(sketchup);
  visual_scenes.appendChild(visual_scene);
  docRootNode.appendChild(visual_scenes);
  docRootNode.appendChild(library_nodes);

  docRootNode.appendChild(library_geometries);

  scene = docNode.createElement('scene');
  instance_visual_scene = docNode.createElement('instance_visual_scene');
  instance_visual_scene.setAttribute('url','#ID2');
  scene.appendChild(instance_visual_scene);
  docRootNode.appendChild(scene);

  f = fopen(filename,'w');
  fprintf(f,'%s',xmlwrite(docNode));
  fclose(f);
end
