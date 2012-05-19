function [Dx,Dy] = readMSH(filename)
  % READMSH Read Adobe Photoshop's .msh format saved from the Liquify tool
  %
  % Inputs:
  %   filename  path to .msh file
  % Outputs:
  %   Dx  h by w x-displacement image 
  %   Dy  h by w y-displacement image 
  %
  fp = fopen(filename,'r');
  % Read HEADER
  header0002 = fread(fp,4,'*uint8');
  if ~all(header0002 == [0;0;0;2])
    warn();
  end
  headerMeshLqfy = fliplr(fread(fp,8,'*char')');
  if ~strcmp(headerMeshLqfy,'MeshLqfy')
    warn();
  end
  header2000 = fread(fp,4,'*uint8');
  if ~all(header2000 == [2;0;0;0])
    warn();
  end
  % Read in size of mesh
  w = fread(fp,1,'*int');
  h = fread(fp,1,'*int');
  assert(w>0);
  assert(h>0);
  % read a zero
  headerzero = fread(fp,1,'*double');
  if headerzero ~= 0
    warn();
  end
  D = fread(fp,2*w*h,'single=>double');
  footerempty = fread(fp,1,'*uint8');
  if ~isempty(footerempty) || ~feof(fp)
    warn();
  end
  fclose(fp);

  Dx = reshape(D(1:2:end),w,h)';
  Dy = reshape(D(2:2:end),w,h)';

  function warn()
    warning('%s might not be a .msh file',filename);
  end

end
