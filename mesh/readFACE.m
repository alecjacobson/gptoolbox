function [F,B] = readFACE(filename)
  % READFACE
  %
  %[F,B] = readFACE(filename)
  %
  % Read triangular faces from a .face file
  % Input:
  %  filename  name of .face file
  % Output:
  %  F  list of triangle indices
  %  B  list of boundary markers
  fp = fopen(filename,'r');
  header = fscanf(fp,'%d %d',2);
  sizeF = header(1);
  boundary_markers = header(2);
  parser = '%d %d %d %d';
  num_items = 4;
  if(boundary_markers ~= 0)
    parser = [parser ' %d'];
    num_items = 5;
  end

  F = fscanf(fp,parser,num_items*sizeF);
  fclose(fp);

  F = reshape(F,num_items,sizeF)';
  B = [];
  if boundary_markers
    B = F(:,5);
  end
  % get rid of indices and boundary markers and make one indexed
  F = F(:,2:4);
end
