function [V, F] = readWRL(filename)
% readWRL - load a mesh from a VRML file
%
%   [V, F] = readWRL(filename);
%
%   Copyright (c) 2004 Gabriel Peyr

fp = fopen(filename,'r');
if fp == -1
 fclose all;
 str = sprintf('Cannot open file %s \n',filename);
 error(str);
end

tempstr = ' ';

key = 'point [';
V = [];

while ( tempstr ~= -1)
  tempstr = fgets(fp);   % -1 if eof
  if( ~isempty(findstr(tempstr,key)) )
    nc = 3;
    while nc>0
      tempstr = fgets(fp);   % -1 if eof
      [cvals,nc] = sscanf(tempstr,'%f %f %f,');
      if nc>0
        if mod(nc,3)~=0
          error('Not correct WRL format');
        end
        V = [V, reshape(cvals, 3, nc/3)];
      end
    end
    break;
  end
end

V=V';

if nargout==1
  return;
end

key = 'coordIndex [';
F = [];

while ( tempstr ~= -1)
  tempstr = fgets(fp);   % -1 if eof
  if( ~isempty(findstr(tempstr,key)) )
    nc = 3;
    while nc>0
      tempstr = fgets(fp);   % -1 if eof
      [cvals,nc] = sscanf(tempstr,'%d %d %d -1,');
      if nc==0 || mod(nc,3)~=0
        [cvals,nc] = sscanf(tempstr,'%d, %d, %d, -1,');
      end
      if nc>0
        F = [F, reshape(cvals, 3, nc/3)+1];
      end
    end
    fclose(fp);
    F = F';
    return;
  end
end

end
