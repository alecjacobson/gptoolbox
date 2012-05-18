function [W] = readDMAT(filename)
  % READDMAT  read a matrix from a dmat file.
  %   first line is <# columns> <# rows>, then values with columns running
  %   faster
  %
  % [W] = readDMAT(filename)
  %
  % Input:
  %   filename  name of .dmat file
  % Output:
  %   W  matrix read from file
  %
  %
  
  % open file
  fp = fopen(filename,'r');
  % read header
  size_W = fliplr(fscanf(fp,'%d %d',[1 2]));
  % read data
  W = fscanf(fp,'%g',size_W);

  if ~feof(fp)
    [size_B,c_B] = fscanf(fp,'%d %d\n',[1 2]);
    size_B = fliplr(size_B);
    if c_B==2
      assert(~any(size_W));
      assert(isempty(W));
      % We're reading binary then
      size_W = size_B;
      [W,cW] = fread(fp,prod(size_W),'*double');
      assert(cW == prod(size_W));
      W = reshape(W,size_W);
    end
  end

  % close file
  fclose(fp);

  % size should match header
  if(~all(size_W == size(W)))
    error('Size in header did not match size of data in file');
  end

end
