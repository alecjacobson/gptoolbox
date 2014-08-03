function [V,F] = readPLY(filename)
  % READPLY Wrapper for read_ply
  % 
  % Input:
  %   filename  path to .ply file
  % Outputs:
  %   V  #V by 3 list of vertex positions
  %   F  #F by 3 list of face indices
  %
  [V,F] = read_ply(filename);

end
