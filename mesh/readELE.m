function E = readELE(ele_file_name)
  % READELE Read (face or tet) elements from a .ele file, used by Stellar and
  % Triangle
  %
  % E = readELE(ele_file_name)
  % 
  % Input:
  %   ele_file_name  path to .ele file
  %
  % Output:
  %   E  list of element indices (1-indexed), could be triangles or tets based
  %     on header of .ele file
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also readNODE
  %

  ele_file_handle = fopen(ele_file_name);
  line = fscanf(ele_file_handle,' %[^\n]s');
  [header,count] = sscanf(line,'%d %d %d',3);
  % number of elements
  num_e = header(1);
  % size of an element
  size_e = header(2);
  % number of attributes
  if(count > 2)
    num_a = header(3);
  else
    num_a = 0;
  end
  if(num_a ~= 0)
    error('Attributes are not supported.');
  end

  % line always start with row index
  parser = '%d';
  num_items = 1;
  if(size_e < 1)
    error(['Size of an element must be at least one, not ' num2str(size_e)]);
  end
  % append to parser enough to read all entries in element
  parser = [parser repmat(' %d',1,size_e)];
  num_items = num_items + size_e;
  E = fscanf(ele_file_handle,parser,[num_items, num_e])';
  fclose(ele_file_handle);

  % get rid of row indices and make one indexed
  E = E(:,2:end) + 1;
end

