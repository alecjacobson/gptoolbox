function writeDAE(filename,V,F)
  % WRITEDAE  Write a mesh to a Collada .dae scene file.
  %
  % writeDAE(filename,V,F)
  %
  % Inputs:
  %   filename  path to .dae file
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of triangle indices
  %

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
  node = docNode.createElement('node');
  node.setAttribute('id','ID3');
  node.setAttribute('name','group_0');

  matrix = docNode.createElement('matrix');
  matrix.appendChild(docNode.createTextNode(sprintf('%d ',eye(4))));
  node.appendChild(matrix);
  instance_geometry = docNode.createElement('instance_geometry');
  instance_geometry.setAttribute('url','#ID4');
  bind_material = docNode.createElement('bind_material');
  bind_material.appendChild(docNode.createElement('technique_common'));
  instance_geometry.appendChild(bind_material);
  node.appendChild(instance_geometry);
  sketchup.appendChild(node);
  visual_scene.appendChild(sketchup);
  visual_scenes.appendChild(visual_scene);
  docRootNode.appendChild(visual_scenes);

  library_geometries = docNode.createElement('library_geometries');
  geometry = docNode.createElement('geometry');
  geometry.setAttribute('id','ID4');
  mesh = docNode.createElement('mesh');
  source = docNode.createElement('source');
  source.setAttribute('id','ID7');
  float_array = docNode.createElement('float_array');
  float_array.setAttribute('id','ID10');
  float_array.setAttribute('count',num2str(numel(V)));
  float_array.appendChild(docNode.createTextNode(sprintf('%g ',V')));
  source.appendChild(float_array);
  technique_common = docNode.createElement('technique_common');
  accessor = docNode.createElement('accessor');
  accessor.setAttribute('count',num2str(size(V,1)));
  accessor.setAttribute('source','#ID10');
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
  vertices.setAttribute('id','ID9');
  input = docNode.createElement('input');
  input.setAttribute('semantic','POSITION');
  input.setAttribute('source','#ID7');
  vertices.appendChild(input);
  mesh.appendChild(vertices);
  triangles = docNode.createElement('triangles');
  triangles.setAttribute('count',num2str(size(F,1)));
  input = docNode.createElement('input');
  input.setAttribute('offset','0');
  input.setAttribute('semantic','VERTEX');
  input.setAttribute('source','#ID9');
  triangles.appendChild(input);
  p = docNode.createElement('p');
  p.appendChild(docNode.createTextNode(sprintf('%d ',(F-1)')));
  triangles.appendChild(p);
  mesh.appendChild(triangles);
  geometry.appendChild(mesh);
  library_geometries.appendChild(geometry);
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
