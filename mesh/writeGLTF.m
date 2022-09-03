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
  %     MT  #V by 3 by #M list of morph targets
  %   'MorphWeights' followed by
  %     MW  #frames by #M list of morph weights
  %
  % Examples:
  %   !python -m json.tool octopus.gltf > octopus-pretty.gltf

  T = [];
  W = [];
  MT = [];
  MW = [];
  N = [];
  fps = 30;
  nf = 0;
  % Map of parameter names to variable names
  params_to_variables = containers.Map( ...
    {'FPS','MorphTargets','MorphWeights','Normals','SkinningTransforms','SkinningWeights'}, ...
    {'fps','MT','MW','N','T','W'});
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


  BYTE = 5120;
  UNSIGNED_BYTE = 5121;
  SHORT = 5122;
  UNSIGNED_SHORT = 5123;
  UNSIGNED_INT = 5125;
  FLOAT = 5126;

  % Number of skinning weights/handles/bones
  m = size(W,2) * ~isempty(W);
  if m>0
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
    %T(:,:,m+(1:k4-k),:) = repmat(eye(3,4),[1 1 k4-k size(T,4)]);
    %m = m+k4-k;

    % Could be unsigned byte or unsigned short. Let's just use byte
    if max(WIk(:)) < 2^8
      WIk_size = 1;
      WIk_type = UNSIGNED_BYTE;
      WIk_matlab_type = 'uint8';
    else
      WIk_size = 2;
      WIk_type = UNSIGNED_SHORT;
      WIk_matlab_type = 'uint16';
    end
    nf = size(T,4) * ~isempty(T);
    if nf>0
      assert(size(T,3) == m);
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
    end
  end
  % Number of morph targets
  mt = size(MT,3) * ~isempty(MT);
  if mt>0
    assert(size(MW,2) == mt);
    nf = size(MW,1);
  end


  if max(F(:)) < 2^16
    F_size = 2;
    F_type = UNSIGNED_SHORT;
    F_matlab_type = 'uint16';
  else
    F_size = 4;
    F_type = UNSIGNED_INT;
    F_matlab_type = 'uint32';
  end


  % Using uint16
  tmpfile = 'writeGLTF.tmp';
  fp = fopen(tmpfile,'wb');
  offsets = 0;

  fwrite(fp,F'-1,F_matlab_type);
  % pad so that next block starts aligned at 4-byte
  Fpad = mod(4-mod(numel(F),4),4);
  fwrite(fp,zeros(Fpad,1),F_matlab_type);
  offsets(end+1) = offsets(end) + (numel(F)+Fpad)*F_size;

  fwrite(fp,V','single');
  % pad so that next block starts aligned at 4-byte
  % Why is this 8 and not 4?
  Vpad = mod(8-mod(numel(V),8),8);
  fwrite(fp,zeros(Vpad,1),'single');
  offsets(end+1) = offsets(end) + (numel(V)+Vpad)*4;

  if ~isempty(N)
    fwrite(fp,N','single');
    % pad so that next block starts aligned at 4-byte
    % Why is this 8 and not 4?
    Npad = mod(8-mod(numel(N),8),8);
    fwrite(fp,zeros(Npad,1),'single');
    offsets(end+1) = offsets(end) + (numel(N)+Npad)*4;
  end

  % Hmm I'm piling this all into a single buffer. Probably it'd be better to put
  % skin weights etc. in their own buffer.
  if m>0
    fwrite(fp,WIk'-1,WIk_matlab_type);
    % pad so that next block starts aligned at 4-byte
    WIkpad = mod(4-mod(numel(WIk),4),4);
    fwrite(fp,zeros(WIkpad,1),WIk_matlab_type);
    offsets(end+1) = offsets(end) + (numel(WIk)+WIkpad)*WIk_size;

    fwrite(fp,Wk','single');
    % pad so that next block starts aligned at 4-byte
    % Why is this 8 and not 4?
    Wkpad = mod(8-mod(numel(Wk),8),8);
    fwrite(fp,zeros(Wkpad,1),'single');
    offsets(end+1) = offsets(end) + (numel(Wk)+Wkpad)*4;
    if nf>0
      % timestamps
      t = (0:nf-1)*1/fps;
      fwrite(fp,t,'single');
      % pad so that next block starts aligned at 4-byte
      % Why is this 8 and not 4?
      tpad = mod(8-mod(numel(t),8),8);
      fwrite(fp,zeros(tpad,1),'single');
      offsets(end+1) = offsets(end) + (numel(t)+tpad)*4;

      fwrite(fp,tT,'single');
      % pad so that next block starts aligned at 4-byte
      % Why is this 8 and not 4?
      tTpad = mod(8-mod(numel(tT),8),8);
      fwrite(fp,zeros(tTpad,1),'single');
      offsets(end+1) = offsets(end) + (numel(tT)+tTpad)*4;

      fwrite(fp,qU,'single');
      % pad so that next block starts aligned at 4-byte
      % Why is this 8 and not 4?
      qUpad = mod(8-mod(numel(qU),8),8);
      fwrite(fp,zeros(qUpad,1),'single');
      offsets(end+1) = offsets(end) + (numel(qU)+qUpad)*4;

      fwrite(fp,sS,'single');
      % pad so that next block starts aligned at 4-byte
      % Why is this 8 and not 4?
      sSpad = mod(8-mod(numel(sS),8),8);
      fwrite(fp,zeros(sSpad,1),'single');
      offsets(end+1) = offsets(end) + (numel(sS)+sSpad)*4;

      fwrite(fp,qVT,'single');
      % pad so that next block starts aligned at 4-byte
      % Why is this 8 and not 4?
      qVTpad = mod(8-mod(numel(qVT),8),8);
      fwrite(fp,zeros(qVTpad,1),'single');
      offsets(end+1) = offsets(end) + (numel(qVT)+qVTpad)*4;
    end
  end
  if mt>0
    for i = 1:mt
      MTi = MT(:,:,i);
      fwrite(fp,MTi','single');
      % pad so that next block starts aligned at 4-byte
      % Why is this 8 and not 4?
      MTipad = mod(8-mod(numel(MTi),8),8);
      fwrite(fp,zeros(MTipad,1),'single');
      offsets(end+1) = offsets(end) + (numel(MTi)+MTipad)*4;
    end
    if nf>0
      % timestamps
      t = (0:nf-1)*1/fps;
      fwrite(fp,t,'single');
      % pad so that next block starts aligned at 4-byte
      % Why is this 8 and not 4?
      tpad = mod(8-mod(numel(t),8),8);
      fwrite(fp,zeros(tpad,1),'single');
      offsets(end+1) = offsets(end) + (numel(t)+tpad)*4;

      fwrite(fp,MW','single');
      % pad so that next block starts aligned at 4-byte
      % Why is this 8 and not 4?
      MWpad = mod(8-mod(numel(MW),8),8);
      fwrite(fp,zeros(MWpad,1),'single');
      offsets(end+1) = offsets(end) + (numel(MW)+MWpad)*4;
    end
  end

  fclose(fp);
  fp = fopen(tmpfile,'r');
  B = fread(fp,offsets(end));
  fclose(fp);
  delete(tmpfile);

  % this expects doubles, but really should be uint8
  Bstr = base64encode(B);
  %assert(isequal(Bstr,matlab.net.base64encode(uint8(B))))

  ARRAY_BUFFER = 34962;
  ELEMENT_ARRAY_BUFFER = 34963;
  bufferViews = {};
  offset_index = 1;
  F_bufferView = numel(bufferViews);
  bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index),'byteLength',numel(F)*F_size,'target',ELEMENT_ARRAY_BUFFER);
  offset_index = offset_index+1;
  V_bufferView = numel(bufferViews);
  bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index),'byteLength',numel(V)*4,'target',ARRAY_BUFFER);
  offset_index = offset_index+1;
  if ~isempty(N)
    N_bufferView = numel(bufferViews);
    bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index),'byteLength',numel(N)*4,'target',ARRAY_BUFFER);
    offset_index = offset_index+1;
  end

  if m>0
    WIk_bufferView = numel(bufferViews);
    bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index+0),'byteLength',numel(WIk)*WIk_size,'target',ARRAY_BUFFER);
    offset_index = offset_index+1;
    Wk_bufferView = numel(bufferViews);
    bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index+1),'byteLength',numel(Wk)*4,'target',ARRAY_BUFFER);
    offset_index = offset_index+1;
    if nf>0
      t_bufferView = numel(bufferViews);
      bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index+0),'byteLength',numel(t)*4);
      tT_bufferView = numel(bufferViews);
      bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index+1),'byteLength',numel(tT)*4);
      qU_bufferView = numel(bufferViews);
      bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index+2),'byteLength',numel(qU)*4);
      sS_bufferView = numel(bufferViews);
      bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index+3),'byteLength',numel(sS)*4);
      qVT_bufferView = numel(bufferViews);
      bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index+4),'byteLength',numel(qVT)*4);
      offset_index = offset_index+5;
    end
  end
  if mt>0
    MT_bufferViews = zeros(mt,1);
    for i = 1:mt
      MT_bufferViews(i) = numel(bufferViews);
      bufferViews{end+1} =  ...
        struct('buffer',0,'byteOffset',offsets(offset_index),'byteLength',numel(MT(:,:,i))*4,'target',ARRAY_BUFFER);
      offset_index = offset_index+1;
    end
    if nf>0
      t_bufferView = numel(bufferViews);
      bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index),'byteLength',numel(t)*4);
      offset_index = offset_index+1;
      MW_bufferView = numel(bufferViews);
      bufferViews{end+1} = struct('buffer',0,'byteOffset',offsets(offset_index),'byteLength',numel(MW)*4);
      offset_index = offset_index+1;
    end
  end

  accessors = {};
  accessors{end+1} = struct('bufferView',F_bufferView,'byteOffset',0,'componentType',F_type,'count',numel(F),'type','SCALAR','max',{{max(F(:))-1}},'min',{{min(F(:))-1}});
  V_accessor = numel(accessors);
  accessors{end+1} = struct('bufferView',V_bufferView,'byteOffset',0,'componentType',FLOAT,'count',size(V,1),'type','VEC3','max',max(V),'min',min(V));
  if ~isempty(N)
    N_accessor = numel(accessors);
    accessors{end+1} = struct('bufferView',N_bufferView,'byteOffset',0,'componentType',FLOAT,'count',size(N,1),'type','VEC3','max',max(N),'min',min(N));
  end
  if m>0
    WIk_accessor = numel(accessors);
    accessors{end+1} = struct('bufferView',WIk_bufferView,'byteOffset',0,'componentType',WIk_type,'count',size(WIk,1),'type','VEC4','max',max(WIk)-1,'min',min(WIk)-1);
    Wk_accessor = numel(accessors);
    accessors{end+1} = struct('bufferView',Wk_bufferView,'byteOffset',0,'componentType',FLOAT,'count',size(Wk,1),'type','VEC4','max',max(Wk),'min',min(Wk));
    if nf>0
      vec = @(X) reshape(X,[],1);
      accessors{end+1} = struct('bufferView',t_bufferView,'byteOffset',0,'componentType',FLOAT,'count',numel(t),'type','SCALAR','max',{{max(t)}},'min',{{min(t)}});
      t_accessor = numel(accessors)-1;
      tT_accessor = [];
      qU_accessor = [];
      sS_accessor = [];
      qVT_accessor = [];
      for i = 1:m
        accessors{end+1} = struct('bufferView',tT_bufferView,'byteOffset',((i-1)*3*nf)*4,'componentType',FLOAT,'count',nf,'type','VEC3','max',vec(max(tT(:,:,:,i),[],3)),'min',vec(min(tT(:,:,:,i),[],3)));
        tT_accessor = [tT_accessor;numel(accessors)-1];
        accessors{end+1} = struct('bufferView',qU_bufferView,'byteOffset',((i-1)*4*nf)*4,'componentType',FLOAT,'count',nf,'type','VEC4','max',vec(max(qU(:,:,i),[],2)),'min',vec(min(qU(:,:,i),[],2)));
        qU_accessor = [qU_accessor;numel(accessors)-1];
        accessors{end+1} = struct('bufferView',sS_bufferView,'byteOffset',((i-1)*3*nf)*4,'componentType',FLOAT,'count',nf,'type','VEC3','max',vec(max(sS(:,:,:,i),[],3)),'min',vec(min(sS(:,:,:,i),[],3)));
        sS_accessor = [sS_accessor;numel(accessors)-1];
        accessors{end+1} = struct('bufferView',qVT_bufferView,'byteOffset',((i-1)*4*nf)*4,'componentType',FLOAT,'count',nf,'type','VEC4','max',vec(max(qVT(:,:,i),[],2)),'min',vec(min(qVT(:,:,i),[],2)));
        qVT_accessor = [qVT_accessor;numel(accessors)-1];
      end
    end
  end
  if mt>0
    MT_accessors = zeros(mt,1);
    for i = 1:mt
      MT_accessors(i) = numel(accessors);
      accessors{end+1} = ...
        struct('bufferView',MT_bufferViews(i),'byteOffset',0,'componentType',FLOAT,'count',size(MT,1),'type','VEC3','max',max(MT(:,:,i)),'min',min(MT(:,:,i)));
    end
    t_accessor = numel(accessors);
    accessors{end+1} = ...
      struct('bufferView',t_bufferView,'byteOffset',0,'componentType',FLOAT,'count',numel(t),'type','SCALAR','max',{{max(t)}},'min',{{min(t)}});
    MW_accessor = numel(accessors);
    accessors{end+1} = ...
      struct('bufferView',MW_bufferView,'byteOffset',0,'componentType',FLOAT,'count',numel(MW),'type','SCALAR','max',{{max(MW(:))}},'min',{{min(MW(:))}});
  end

  buffer = struct( ...
    'uri',['data:application/octet-stream;base64,' Bstr], ...
    'byteLength', numel(B));

  % First node is the mesh
  attributes = struct("POSITION",V_accessor);
  if ~isempty(N)
    attributes.NORMAL = N_accessor;
  end
  if m>0
    attributes.JOINTS_0 = WIk_accessor;
    attributes.WEIGHTS_0 = Wk_accessor;
  end
  mesh = struct( ...
    "primitives",{{struct("attributes",attributes,"indices",0)}});
  if mt>0
    mesh.primitives{1}.targets = {};
    for i = 1:mt
      mesh.primitives{1}.targets{i} = struct('POSITION',MT_accessors(i));
    end
    mesh.weights = zeros(1,mt);
  end
  nodes_list = {0};
  nodes = {struct("mesh",0)};

  if m>0
    nodes{end}.skin = 0;
    %% Put all the nodes at the root level
    %nodes_list(end+(1:m)) = num2cell(nodes_list{end}+(1:m));


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
  if m>0
    gltf.skins = {struct('joints',{num2cell(0+1:m)})};
    if nf>0
      samplers = {};
      channels = {};
      for i = 1:m
        samplers{end+1} = struct('input',t_accessor,'output',tT_accessor(i),'interpolation', 'STEP');
        samplers{end+1} = struct('input',t_accessor,'output',qU_accessor(i),'interpolation', 'STEP');
        samplers{end+1} = struct('input',t_accessor,'output',sS_accessor(i),'interpolation', 'STEP');
        samplers{end+1} = struct('input',t_accessor,'output',qVT_accessor(i),'interpolation','STEP');

        channels{end+1} = struct('sampler',4*(i-1)+0,'target',struct('node',m+i,'path','translation'));
        channels{end+1} = struct('sampler',4*(i-1)+1,'target',struct('node',m+i,'path','rotation'));
        channels{end+1} = struct('sampler',4*(i-1)+2,'target',struct('node',m+i,'path','scale'));
        channels{end+1} = struct('sampler',4*(i-1)+3,'target',struct('node',i,'path','rotation'));
      end
      gltf.animations = {struct('samplers',{samplers},'channels',{channels})};
    end
  end
  if mt>0
    if nf>0
      samplers = {struct('input',t_accessor,'output',MW_accessor,'interpolation','LINEAR')};
      channels = {struct('sampler',0,'target',struct('node',0,'path','weights'))};
      gltf.animations = {struct('samplers',{samplers},'channels',{channels})};
    end
  end
  fid = fopen(filename,'w');
  fprintf(fid,'%s',jsonencode(gltf));
  fclose(fid);

  %tic;
  %[s,r] = system(sprintf('cat %s | ruby -r json -e "jj JSON.parse gets"',filename));r
  %toc



end
