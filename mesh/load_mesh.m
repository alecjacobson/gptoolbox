function [V,F] = load_mesh(filename)
  % read in vertices and faces from a .off or .obj file
  % Input:
  %   filename  file holding mesh
  % Output:
  %   V (vertex list) 
  %   F (face list) fields
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: readOBJ, readOBJfast, readOFF
  %
  if ~isempty(regexp(filename,'\.off$'))
    [V,F] = readOFF(filename);
  elseif ~isempty(regexp(filename,'\.obj$'))
    try
      [V,F] = readOBJfast(filename);
    catch exception
      fprintf('Fast reader failed, retrying with more robust, slower reader\n');
      [V,F] = readOBJ(filename);
    end
  else
    error('Input file must be .off or .obj file.');
  end
end
