function [V,E,P,BE,CE,PE] = readTGF(filename)
  % READTGF
  %
  % [V,E,P,BE,CE,PE] = readTGF(filename)
  %
  % Read a graph from a .tgf file
  %
  % Input:
  %  filename  .tgf file name
  % Ouput:
  %  V  # vertices by 3 list of vertex positions
  %  E  # edges by 2 list of edge indices
  %  P  # point-handles list of point handle indices
  %  BE # bone-edges by 2 list of bone-edge indices
  %  CE # cage-edges by 2 list of cage-edge indices
  %  PE # pseudo-edges by 2 list of pseudo-edge indices
  % 
  % Assumes that graph vertices are 3 dimensional

  V = [];
  E = [];
  BE = [];
  CE = [];
  PE = [];
  fp = fopen(filename,'r');
  % read next whole line
  line = fscanf(fp,' %[^\n]s');
  % read vertices until seeing separator
  while(line(1,1) ~= '#')
    [v,count] = sscanf(line,'%d %g %g %g',4);
    % first number is always index
    v = v(2:min(count,4));
    V = [V;v'];
    line = fscanf(fp,' %[^\n]s');
  end
  % read edges until file is through
  line = fscanf(fp,' %[^\n]s');
  while(sum(size(line)) > 0 && line(1,1) ~= '#')
    [e,count] = sscanf(line,'%d %d %d %d %d',5);
    assert(count >= 2);
    % add edge
    % .tgf format is 1-indexed
    E = [E;e(1:2)'];
    if count >= 3 && e(3)
      BE = [BE;e(1:2)'];
    end
    % bone edges trump pseudo-edges
    if count >= 4 && e(4) && ~e(3)
      PE = [PE;e(1:2)'];
    end
    % bone edges trump cage edges
    if count >= 5 && e(5) && ~e(3)
      CE = [CE;e(1:2)'];
    end
    line = fscanf(fp,' %[^\n]s');
  end

  % point handles are points not attached to a bone
  P = 1:size(V,1);
  P = P(~ismember(P,BE(:)));
  fclose(fp);
end
