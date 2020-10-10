function [V,WI,P] = readBF(filename)
  % READBF Read a skeleton from the .bf bone forest format
  %
  % [V,WI,P] = readBF(filename)
  %
  % Read bone forest from a .bf file
  %
  % Input:
  %  filename  .bf file name
  % Output:
  %  V  # vertices by 3 list of vertex offsets from parents
  %  WI  # vertices list of indices of weights (0 means no weights)
  %  P  # vertices list of indices of bone parents (0 means root)
  % 
  D = dlmread(filename);
  V = D(:,3:end);
  WI = D(:,1)+1;
  P = D(:,2)+1;
end
