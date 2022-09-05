function writeGLTF(filename,V,F,varargin)
  % writeGLTF(filename,V,F)
  %
  % Inputs:
  %   filename  path to .gltf file
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices into rows of V
  %   Optional:
  %   'SkinningTransforms' followed by
  %     T  3 by 4 by #T by #A list of animated global skinning transformations
  %   'SkinningWeights' followed by
  %     W  #V by #T list of skinning weights
  %   'MorphTargets' followed by
  %     MV  #V by 3 by #M list of morph targets
  %   'MorphWeights' followed by
  %     MW  #frames by #M list of morph weights
  %   'TextureCoordinates' followed by
  %     UV  #V by 2 list of texture coordinates
  %
  % Examples:
  %   !python -m json.tool octopus.gltf > octopus-pretty.gltf
  NEAREST = 9728;
  LINEAR = 9729;
  NEAREST_MIPMAP_NEAREST = 9984;
  LINEAR_MIPMAP_NEAREST = 9985;
  NEAREST_MIPMAP_LINEAR = 9986;
  LINEAR_MIPMAP_LINEAR = 9987;
  CLAMP_TO_EDGE = 33071;
  MIRRORED_REPEAT = 33648;
  REPEAT = 10497;

  T = [];
  W = [];
  MV = [];
  MN = [];
  MW = [];
  N = [];
  fps = 30;
  nf = 0;
  UV = [];
  tex_im = [];
  mag_filter = LINEAR;
  min_filter = LINEAR_MIPMAP_LINEAR;
  wrap_s = MIRRORED_REPEAT;
  wrap_t = MIRRORED_REPEAT;


  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'FPS','MagFilter','MinFilter','MorphTargets','MorphNormals','MorphWeights','Normals','SkinningTransforms','SkinningWeights','TextureCoordinates','TextureImage'}, ...
    {'fps','mag_filter','min_filter','MV','MN','MW','N','T','W','UV','tex_im'});
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

  mag_filter = string_to_magic_value(mag_filter);
  min_filter = string_to_magic_value(min_filter);
  function value = string_to_magic_value(string)
    if ischar(string)
      switch(string)
      case 'LINEAR'
        value = LINEAR;
      case 'NEAREST'
        value = NEAREST;
      case 'LINEAR_MIPMAP_NEAREST'
        value = LINEAR_MIPMAP_NEAREST;
      case 'LINEAR_MIPMAP_LINEAR'
        value = LINEAR_MIPMAP_LINEAR;
      otherwise
        error([string ' not found/implemented.']);
      end
    else
      value = string;
    end
  end


  % set flags 
  has_normals = ~isempty(N);
  has_skinning_weights = ~isempty(W);
  has_skinning_transforms = ~isempty(T);
  has_morph_target_positions = ~isempty(MV);
  has_morph_target_weights = ~isempty(MW);
  has_morph_target_normals = ~isempty(MN);
  has_texture_coordinates = ~isempty(UV);
  has_texture_image = ~isempty(tex_im);



  % Determine what type to use for faces
  BYTE = 5120;
  UNSIGNED_BYTE = 5121;
  SHORT = 5122;
  UNSIGNED_SHORT = 5123;
  UNSIGNED_INT = 5125;
  FLOAT = 5126;
  ARRAY_BUFFER = 34962;
  ELEMENT_ARRAY_BUFFER = 34963;
  if max(F(:)) < 2^16
    F_size = 2;
    F_type = UNSIGNED_SHORT;
    F_matlab_type = 'uint16';
  else
    F_size = 4;
    F_type = UNSIGNED_INT;
    F_matlab_type = 'uint32';
  end

  % prepare things to write
  buffered_data = {};
  buffered_data{end+1} = struct('Name','F','Data',F,'Type','SCALAR','ComponentType',F_type,'Size',F_size,'MatlabType',F_matlab_type,'Target',ELEMENT_ARRAY_BUFFER);
  buffered_data{end+1} = struct('Name','V','Data',V,'Target',ARRAY_BUFFER);
  if has_normals
    buffered_data{end+1} = struct('Name','N','Data',normalizerow(N),'Target',ARRAY_BUFFER);
  end
  if has_morph_target_weights
    assert(size(MW,2) == size(MV,3));
    buffered_data{end+1} = struct('Name','MW','Data',MW,'Type','SCALAR');
    MWt = (0:size(MW,1)-1)*1/fps;
    buffered_data{end+1} = struct('Name','MWt','Data',MWt,'Type','SCALAR');
  end
  if has_morph_target_positions
    buffered_data{end+1} = struct('Name','MV','Data',MV,'Target',ARRAY_BUFFER);
  end
  if has_morph_target_normals
    buffered_data{end+1} = struct('Name','MN','Data',MN,'Target',ARRAY_BUFFER);
  end
  if has_skinning_weights
    [WIk,Wk] = prepare_skinning_weights(W);
    if max(WIk(:)) < 2^8
      WIk_size = 1;
      WIk_type = UNSIGNED_BYTE;
      WIk_matlab_type = 'uint8';
    else
      WIk_size = 2;
      WIk_type = UNSIGNED_SHORT;
      WIk_matlab_type = 'uint16';
    end
    buffered_data{end+1} = struct('Name','WIk','Data',WIk,'Size',WIk_size,'MatlabType',WIk_matlab_type,'ComponentType',WIk_type,'Target',ARRAY_BUFFER);
    buffered_data{end+1} = struct('Name','Wk','Data',Wk,'Target',ARRAY_BUFFER);
  end
  if has_skinning_transforms
    assert(size(W,2) == size(T,3));
    [tT,qU,sS,qVT] = prepare_skinning_transforms(T);
    buffered_data{end+1} = struct('Name','tT','Data',tT);
    buffered_data{end+1} = struct('Name','qU','Data',qU);
    buffered_data{end+1} = struct('Name','sS','Data',sS);
    buffered_data{end+1} = struct('Name','qVT','Data',qVT);
    Tt = (0:size(T,4)-1)*1/fps;
    buffered_data{end+1} = struct('Name','Tt','Data',Tt,'Type','SCALAR');
  end
  if has_texture_coordinates
    buffered_data{end+1} = struct('Name','UV','Data',UV,'Target',ARRAY_BUFFER);
  end

  [buffer,bufferViews,accessors,accessor_hash] = prepare_buffer(buffered_data);

  % First node is the mesh
  attributes = struct("POSITION",accessor_hash.V);

  if has_normals
    attributes.NORMAL = accessor_hash.N;
  end
  if has_texture_coordinates
    attributes.TEXCOORD_0 = accessor_hash.UV;
  end
  if has_skinning_weights
    attributes.JOINTS_0 = accessor_hash.WIk;
    attributes.WEIGHTS_0 = accessor_hash.Wk;
  end
  mesh = struct( ...
    "primitives",{{struct("attributes",attributes,"indices",accessor_hash.F)}});
  if has_morph_target_positions
    if ~isfield(mesh.primitives{1},'targets')
      mesh.primitives{1}.targets = num2cell(repmat(struct(),size(MV,3),1));
    end
    for i = 1:size(MV,3)
      mesh.primitives{1}.targets{i} = setfield(mesh.primitives{1}.targets{i},'POSITION',accessor_hash.MV(i));
    end
    mesh.weights = zeros(1,size(MV,3));
  end
  if has_morph_target_normals
    if ~isfield(mesh.primitives{1},'targets')
      mesh.primitives{1}.targets = num2cell(repmat(struct(),size(MN,3),1));
    end
    for i = 1:size(MN,3)
      mesh.primitives{1}.targets{i} = setfield(mesh.primitives{1}.targets{i},'NORMAL',accessor_hash.MN(i));
    end
  end

  nodes_list = {0};
  nodes = {struct("mesh",0)};
  if has_skinning_weights
    m = size(W,2);
    nodes{end}.skin = 0;
    % First put leaf nodes that will animate V-rotation
    for i = 1:m
      nodes{end+1} = struct( ...
        "name",sprintf('V-%d',i), ...
        "rotation",[0 0 0 1]);
    end
    % Then put parent auxiliary nodes that will animate translation, U-rotation, scale
    for i = 1:m
      nodes{end+1} = struct( ...
        "name",sprintf('TUS-%d',i), ...
        "children",{{i}}, ...
        "translation",[0 0 0], ...
        "rotation",[0 0 0 1], ...
        "scale",[1 1 1]);
    end
    %Put fake root at scene level and attach parents
    nodes_list{end+1} = 2*m+1;
    common_root = struct("name","root","children",{num2cell(m+(1:m))});
    nodes = {nodes{:} common_root};
  end

  % These {{double-brackets}} are needed when using the `struct` function.
  gltf = struct( ...
    "scene",0, ...
    "scenes",{{struct("nodes",{nodes_list})}}, ...
    "nodes",{nodes}, ...
    "meshes",{{mesh}}, ...
    "buffers",{{buffer}}, ...
    "bufferViews",{bufferViews}, ...
    "accessors",{accessors}, ...
    "asset",struct("version","2.0") ...
    );

  % Are these just for animation or do they overlap with texture samplers?
  samplers = {};
  channels = {};
  if has_skinning_weights
    gltf.skins = {struct('joints',{num2cell(0+1:m)})};
    if has_skinning_transforms
      for i = 1:m
        channels{end+1} = struct('sampler',numel(samplers),'target',struct('node',m+i,'path','translation'));
        samplers{end+1} = struct('input',accessor_hash.Tt,'output',accessor_hash.tT(i),'interpolation', 'STEP');
        channels{end+1} = struct('sampler',numel(samplers),'target',struct('node',m+i,'path','rotation'));
        samplers{end+1} = struct('input',accessor_hash.Tt,'output',accessor_hash.qU(i),'interpolation', 'STEP');
        channels{end+1} = struct('sampler',numel(samplers),'target',struct('node',m+i,'path','scale'));
        samplers{end+1} = struct('input',accessor_hash.Tt,'output',accessor_hash.sS(i),'interpolation', 'STEP');
        channels{end+1} = struct('sampler',numel(samplers),'target',struct('node',i,'path','rotation'));
        samplers{end+1} = struct('input',accessor_hash.Tt,'output',accessor_hash.qVT(i),'interpolation','STEP');
      end
    end
  end
  if has_morph_target_weights
    channels{end+1} = struct('sampler',numel(samplers),'target',struct('node',0,'path','weights'));
    samplers{end+1} = struct('input',accessor_hash.MWt,'output',accessor_hash.MW,'interpolation','LINEAR');
  end
  if ~isempty(samplers)
    gltf.animations = {struct('samplers',{samplers},'channels',{channels})};
  end
  if has_texture_image 

    if mag_filter == LINEAR
      tmpfile = 'writeGLTF.jpg';
      imwrite(flipud(tex_im),tmpfile,'Quality',100);
    else
      tmpfile = 'writeGLTF.png';
      imwrite(flipud(tex_im),tmpfile);
    end
    tex_image_uri = imdata(tmpfile);
    delete(tmpfile);


    gltf.images = {struct('uri',tex_image_uri)};
    tex_image_index = numel(gltf.images)-1;
    if ~isfield(gltf,'samplers')
      gltf.samplers = {};
    end
    gltf.samplers{end+1} = struct('magFilter',mag_filter,'minFilter',min_filter,'wrapS',wrap_s,'wrapT',wrap_t);
    gltf.textures = {struct('sampler',numel(gltf.samplers)-1,'source',tex_image_index)};

    if has_texture_coordinates
      gltf.materials = {struct( ...
        'name','default', ...
        'pbrMetallicRoughness', ...
          struct( ...
            'baseColorTexture',struct('index',0), ...
            'metallicFactor',0, ...
            'roughnessFactor',1))};
      gltf.meshes{1}.primitives{1}.material = numel(gltf.materials)-1;
    end
  end

  fid = fopen(filename,'w');
  fprintf(fid,'%s',jsonencode(gltf));
  fclose(fid);

  function [buffer,bufferViews,accessors,accessor_hash] = prepare_buffer(buffered_data)
    buffer_index = 0;
    bufferViews = {};
    accessors = {};
    accessor_hash = struct();
    byte_count = 0;
    % write all data, tracking offsets
    tmpfile = 'writeGLTF.tmp';
    fp = fopen(tmpfile,'wb');
    for bi = 1:numel(buffered_data)
      item = buffered_data{bi};
      if isfield(item,'Size')
        adjusted_data = item.Data - 1;
        pad_amount = 4;
      else
        % default to single
        item.ComponentType = FLOAT;
        item.Size = 4;
        item.MatlabType = 'single';
        pad_amount = 8;
        adjusted_data = item.Data;
      end
      if ~isfield(item,'Type')
        switch(size(adjusted_data,2))
        case 2
          item.Type = 'VEC2';
        case 3
          item.Type = 'VEC3';
        case 4
          item.Type = 'VEC4';
        otherwise
          item.Type = 'SCALAR';
        end
      end

      accessor_ids = numel(accessors)-1 + (1:size(item.Data,3));
      for i = 1:size(item.Data,3)
        adjusted_data_i = adjusted_data(:,:,i);
        switch(item.Type)
        case 'SCALAR'
          accessor_count = numel(adjusted_data_i);
          accessor_min = {{min(  adjusted_data_i,[],[1 2])}};
          accessor_max = {{max(  adjusted_data_i,[],[1 2])}};
        case {'VEC2','VEC3','VEC4'}
          accessor_count = size( adjusted_data_i,1);
          accessor_min = min(    adjusted_data_i,[],1);
          accessor_max = max(    adjusted_data_i,[],1);
        end

        fwrite(fp,adjusted_data_i',item.MatlabType);
        % pad so that next block starts aligned at 4-byte
        bi_pad = mod(pad_amount-mod(numel(adjusted_data_i),pad_amount),pad_amount);
        fwrite(fp,zeros(bi_pad,1),item.MatlabType);
        % prepare bufferView and accessor
        % Always use 'buffer',0  → caller MUST adjust all if needed
        bufferViews{end+1} = struct( ...
          'name',item.Name, ...
          'buffer',buffer_index, ...
          'byteOffset',byte_count, ...
          'byteLength',numel(adjusted_data_i)*item.Size);
        if isfield(item,'Target')
          bufferViews{end} = setfield(bufferViews{end},'target',item.Target);
        end
        accessors{end+1} = struct( ...
          'name',item.Name, ...
          'bufferView',numel(bufferViews)-1, ...
          'byteOffset',0, ...
          'componentType',item.ComponentType, ...
          'type',item.Type, ...
          'count',accessor_count, ...
          'min',accessor_min, ...
          'max',accessor_max);
        % Increment byte_count
        byte_count = byte_count + (numel(adjusted_data_i)+bi_pad)*item.Size;
      end
      accessor_hash = setfield(accessor_hash,item.Name,accessor_ids);
    end
    fclose(fp);
    fp = fopen(tmpfile,'r');
    B = fread(fp,byte_count);
    fclose(fp);
    delete(tmpfile);
    Bstr = base64encode(B);
    %assert(isequal(Bstr,matlab.net.base64encode(uint8(B))))
    buffer = struct( ...
      'uri',['data:application/octet-stream;base64,' Bstr], ...
      'byteLength', numel(B));
  end

  function [WIk,Wk] = prepare_skinning_weights(W)
    % Number of skinning weights/handles/bones
    m = size(W,2) * ~isempty(W);
    % Convert weight matrix into multiple of 4 index-value pairs
    assert(size(W,2) == m);
    k = max(sum(W~=0,2));
    k4 = ceil(k/4)*4;
    % I'm not sure how >4 is handled.
    assert(k4 == 4);
    [Wk,WIk] = maxk(W,k,2);
    % Pad with null bones
    Wk(:,k+1:k4)  = 0;
    WIk(:,k+1:k4) = 0+1;%repmat(m+(1:k4-k),size(WIk,1),1);
    WIk(Wk==0) = 0+1;
  end

  function [tT,qU,sS,qVT] = prepare_skinning_transforms(T)
    m = size(T,3);
    nf = size(T,4) * ~isempty(T);
    T = permute(T,[1 2 4 3]);
    assert(size(T,2) == 4);
    tT = T(:,4,:,:);
    [sU,sS,sV] = pagesvd(T(1:3,1:3,:,:));
    % Push all of the reflection into sS
    % matlabFunction(det(sym('m',[3 3])))
    pagedet = @(m1_1,m1_2,m1_3,m2_1,m2_2,m2_3,m3_1,m3_2,m3_3) m1_1.*m2_2.*m3_3-m1_1.*m2_3.*m3_2-m1_2.*m2_1.*m3_3+m1_2.*m2_3.*m3_1+m1_3.*m2_1.*m3_2-m1_3.*m2_2.*m3_1;
    pagedet = @(m) pagedet(m(1,1,:,:),m(1,2,:,:),m(1,3,:,:),m(2,1,:,:),m(2,2,:,:),m(2,3,:,:),m(3,1,:,:),m(3,2,:,:),m(3,3,:,:));
    sUdet = sign(pagedet(sU));
    sVdet = sign(pagedet(sV));
    sU = sU.*sUdet;
    sV = sV.*sVdet;
    sS = sS.*sVdet.*sUdet;
    sS = cat(2,sS(1,1,:,:),sS(2,2,:,:),sS(3,3,:,:));

    qU = reshape(mat2quat(reshape(sU,[3 3 m*nf]))',[4 nf m]);
    qVT = reshape(mat2quat(reshape(permute(sV,[2 1 3 4]),[3 3 m*nf]))',[4 nf m]);
    for i = 1:m
      for f = 1:nf
        assert(norm(T(:,:,f,i) - [quat2mat(qU(:,f,i)')*diag(sS(:,:,f,i))*quat2mat(qVT(:,f,i)') tT(:,:,f,i)],inf)<1e-12);
      end
    end
    % gltf uses ijkw
    qU = qU([2 3 4 1],:,:);
    qVT = qVT([2 3 4 1],:,:);

    tT = permute(tT,[3 1 4 2]);
    sS = permute(sS,[3 2 4 1]);
    qU = permute(qU,[2 1 3]);
    qVT = permute(qVT,[2 1 3]);
  end

end
