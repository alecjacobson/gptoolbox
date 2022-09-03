function [V,F] = readOBJ_sequence(prefix)
  % [V,F] = readOBJ_sequence(prefix)
  % 
  % Inputs:
  %   prefix  glob style paths (e.g., 'mymesh.00*.obj')
  % Outputs:
  %   V  #files cell array of mesh vertex positions
  %   F  #files cell array of mesh faces indices into corresponding vertex
  %     positions
  %   ** Whenever possible, cell arrays are reduced into arrays: i.e., if all
  %   elements of V have the same matrix size, then a #V by 3 by #files matrix
  %   is output instead (similarly for F).
  %   ** Whenever possible, duplicated information across files is removed:
  %   i.e., if all elements of F are the same for each file, then a single #F by
  %   3 matrix is returned.
  %
  [status,res] = system(sprintf('printf "%%s\n" %s',prefix));
  files = strsplit(res,'\n');
  files = files(~cellfun(@isempty,files));
  V = cell(1,numel(files));
  F = cell(1,numel(files));
  parfor fi = 1:numel(files)
    progressbar(fi,numel(files));
    [V{fi},F{fi}] = load_mesh(files{fi});
  end
  V = convert_to_array_if_possible(V);
  F = convert_to_array_if_possible(F);
  function V = convert_to_array_if_possible(V)
    if numel(V) > 0 && ...
      all(cell2mat(cellfun(@size,V','UniformOutput',false))==size(V{1}),[1 2])
      V = cell2mat(permute(V,[1 3 2]));
      if size(V,3) > 1 &&  ...
        V(1,1,1)==V(1,1,2) && ...
        all(V(:,:,1)==V(:,:,2),[1 2]) && ...
        all(V(:,:,1)==V(:,:,3:end),[1 2 3])
        V = V(:,:,1);
      end
    end
  end
end
