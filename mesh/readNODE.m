function [V,I] = readNODE(filename)
  % READNODE Read vertex positions from .node file, .node files are used by
  % Stellar and Triangle
  %
  % V = readNODE(filename)
  %
  % Input:
  %  filename  path to .node file
  %
  % Output:
  %  V  list of vertex positions
  %  I  list of indices (first tells whether 0 or 1 indexed)
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: readELE
  %

  node_file_handle = fopen(filename);
  % read header which contains number of vertices, number of dimensions, number
  % of attributes, and number of boundary markers
  header = fscanf(node_file_handle,'%d %d %d %d',4);
  num_v = header(1);
  dim = header(2);
  num_a = header(3);
  num_b = header(4);
  % assemble fscanf string
  % line always contains index 
  parser = '%d ';
  num_items = 1;
  % next each line contains coordinates (as many as there are dimensions)
  if(dim == 2)
    parser = [parser '%lf %lf '];
    num_items = num_items + 2;
  elseif(dim == 3)
    parser = [parser '%lf %lf %lf '];
    num_items = num_items + 3;
  else
    error(['Dimension must be 2 or 3, not ' num2str(dim)]);
  end
  if(num_a ~= 0)
    error('Attributes are not supported.');
  end
  % if there are going to be boundary markers then add spot in parser
  if(num_b == 1)
    parser = [parser '%d'];
    num_items = num_items + 1;
  elseif(num_b ~= 0)
    error(['Number of boundary markers must be 1 or 0, not ' num2str(num_b)]);
  end
  % read num_v lines
  V = fscanf(node_file_handle,parser,[num_items, num_v])';
  fclose(node_file_handle);
  I = V(:,1);
  % lose indices 
  V = V(:,2:end);
  % lose boundary markers
  if(num_b == 1)
    V = V(:,1:(end-1));
  end
end

