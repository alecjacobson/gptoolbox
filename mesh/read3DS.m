function [V,F,flag,TC] = read3DS(filename)
  % READ3DS Read a triangle mesh from a .3ds file. 
  % 
  % [V,F,flag] = read3DS(filename)
  %
  % Inputs:
  %   filename  path to .3ds file
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices into V
  %   flag  #F list of face "flags"
  %   TC  #TC by 2 list of texture coordinates
  %
  % See also: load_mesh
  %

  % http://www.gamedev.net/topic/313126-3ds-parsing-tutorial/
  f = fopen(filename, 'r');

  vertices_id = 16656;
  texture_coords_id = 16704;
  faces_id = 16672;
  mesh_id = 16384;
  traverse = { ...
    19789 ... main
    15677 ... 3d editor
    mesh_id ... mesh
    16640 ... triangle mesh
    vertices_id ...
    texture_coords_id ...
    faces_id ...
    };
  skip = { ...
    2 ... version 
    15678 ... mesh version
    15678 ... edit material
    };

  V = [];
  F = [];
  TC = [];
  depth = 0;
  counts = [];
  lens = [];
  while feof(f)~=1
    id = fread(f,1,'ushort');
    id = uint16(id);
    if isempty(id) || id == 0
      break;
    end
    len = fread(f,1,'int');
    did_traverse = false;
    switch id
    case traverse
      % do something 
      depth = depth + 1;
      lens = [lens len];
      counts = [counts 6];
      did_traverse = true;
      %fprintf('%sTraversing chunk with id %d (%04x) and length %d\n',repmat(' ',1,depth),id,swapbytes(id),len);
      switch id
      case mesh_id
        name = [];
        while true
          ch = fread(f,1,'*char')';
          if ch == 0
            break;
          end
          name = [name ch];
        end
      case vertices_id
        n = fread(f,1,'ushort');
        v = fread(f,n*3,'float32');
        last_V_len = size(V,1);
        V = [V;reshape(v,3,n)'];
        counts(end) = counts(end)+n*3*4+2;
        scatter3(V(:,1),V(:,2),V(:,3),'.');
        axis equal;
      case texture_coords_id
        n = fread(f,1,'ushort');
        tc = fread(f,n*2,'float32');
        counts(end) = counts(end)+n*2*4+2;
        TC = [TC;reshape(tc,2,n)'];
      case faces_id
        n = fread(f,1,'ushort');
        faces = fread(f,n*4,'ushort');
        counts(end) = counts(end)+n*4*2+2;
        F = [F;last_V_len+reshape(faces,4,n)'+1];
      otherwise
        %?
      end
    otherwise
      %fprintf('%sSkipping chunk with id %d (%x) and length %d\n',repmat(' ',1,depth),id,id,len);
      fseek(f,len-6,'cof');
      counts(end) = counts(end)+len;
    end
    %fprintf('%d ',lens);
    %fprintf('\n');
    %fprintf('%d ',counts);
    %fprintf('\n');
    while counts(end) == lens(end)
      depth = depth -1;
      lens = lens(1:end-1);
      counts = [counts(1:end-2) counts(end-1)+counts(end)];
      %fprintf('%s.\n',repmat(' ',1,depth));
    end
  end
  flag = F(:,4);
  F = F(:,1:3);
end
