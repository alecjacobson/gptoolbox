function writeGLTF(filename,V,F)
  % writeGLTF(filename,V,F)
  %
  % Inputs:
  %   filename  path to .gltf file
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh triangle indices into rows of V
  %


  BYTE = 5120;
  UNSIGNED_BYTE = 5121;
  SHORT = 5122;
  UNSIGNED_SHORT = 5123;
  UNSIGNED_INT = 5125;
  FLOAT = 5126;

  if max(F(:)) < 2^16
    F_size = 2;
    F_type = UNSIGNED_SHORT;
    F_matlab_type = 'uint16';
  else
    F_size = 4;
    F_type = UNSIGNED_INT;
    F_matlab_type = 'uint32';
  end

  mesh = struct( ...
    "primitives",{{struct("attributes",struct("POSITION",1),"indices",0)}});



  % Using uint16
  tmpfile = 'writeGLTF.tmp';
  fp = fopen(tmpfile,'wb');
  fwrite(fp,F'-1,F_matlab_type);
  % pad so that next block starts aligned at 4-byte
  Fpad = mod(4-mod(numel(F),4),4);
  fwrite(fp,zeros(Fpad,1),F_matlab_type);
  fwrite(fp,V','single');
  % pad so that next block starts aligned at 4-byte
  Vpad = mod(8-mod(numel(V),8),8);
  fwrite(fp,zeros(Vpad,1),'single');
  fclose(fp);
  fp = fopen(tmpfile,'r');
  B = fread(fp,(numel(F)+Fpad)*F_size+numel(V)*4);
  fclose(fp);
  delete(tmpfile);

  Bstr = base64encode(B);
  %assert(isequal(Bstr,matlab.net.base64encode(B)))


  ARRAY_BUFFER = 34962;
  ELEMENT_ARRAY_BUFFER = 34963;
  bufferViews = [ ...
    struct('buffer',0,'byteOffset',0,'byteLength',numel(F)*F_size,'target',ELEMENT_ARRAY_BUFFER) ...
    struct('buffer',0,'byteOffset',(numel(F)+Fpad)*F_size,'byteLength',numel(V)*4,'target',ARRAY_BUFFER) ...
    ];

  accessors = [ ...
    struct('bufferView',0,'byteOffset',0,'componentType',F_type,'count',numel(F),'type','SCALAR','max',{{max(F(:))-1}},'min',{{min(F(:))-1}}) ...
    struct('bufferView',1,'byteOffset',0,'componentType',FLOAT,'count',size(V,1),'type','VEC3','max',max(V),'min',min(V)) ...
    ];

  buffer = struct( ...
    'uri',['data:application/octet-stream;base64,' Bstr], ...
    'byteLength', numel(B));

  % These {{double-brackets}} are needed when using the `struct` function.
  gltf = struct( ...
    "scene",0, ...
    "scenes",{{struct("nodes",{{0}})}}, ...
    "nodes",{{struct("mesh",0)}}, ...
    "meshes",{{mesh}}, ...
    "buffers",{{buffer}}, ...
    "bufferViews",bufferViews, ...
    "accessors",accessors, ...
    "asset",struct("version","2.0") ...
    );
  fid = fopen(filename,'w');
  fprintf(fid,'%s',jsonencode(gltf));
  fclose(fid);

  %tic;
  %[s,r] = system(sprintf('cat %s | ruby -r json -e "jj JSON.parse gets"',filename));r
  %toc





end
