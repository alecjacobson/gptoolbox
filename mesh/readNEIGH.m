function N = readNEIGH(filename)
  % READNEIGH  Read tetrahedral neighbor information from a .neigh file (as
  % produced by tetgen)
  % 
  % N = readNEIGH(filename)
  %
  % Inputs:
  %   filename  path to .neigh file
  % Outputs:
  %   N  #simplices by #size-of-simplex neighborhood information (-1) indicates
  %     boundary. T(i,j) *should* indicate the neighbor to the jth face of the
  %     ith tet. *However* tetgen does not seem consistent. Consider
  %     post-processing with fixNEIGH.m
  %
  % See also: tt, tetgen, readNODE, readELE, fixNEIGH
  %


  fp = fopen(filename);
  line = fscanf(fp,' %[^\n]s');
  [header,count] = sscanf(line,'%d %d',2);
  if count~=2
    fclose(fp);
    error('Bad header');
  end

  % number of elements
  n = header(1);
  % size of an element
  size_e = header(2);

  parser = '%d';
  % append to parser enough to read all entries in element + 1 for index
  parser = [parser repmat(' %d',1,size_e+1)];
  N = fscanf(fp,parser,[size_e+1 n])';
  fclose(fp);

  % get rid of row indices and make one indexed
  N = N(:,2:end) + 1;
  N(~N) = -1;
end


