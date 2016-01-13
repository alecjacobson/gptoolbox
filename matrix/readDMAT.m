function [W] = readDMAT(filename)
  % READDMAT  read a matrix from a dmat file.  first line is <# columns> <#
  % rows>, then values with columns running faster
  %
  % [W] = readDMAT(filename)
  %
  % Input:
  %   filename  name of .dmat file
  % Output:
  %   W  matrix read from file
  %
  % See also: writeDMAT
  %
  
  % open file
  fp = fopen(filename,'r');
  % read header
  size_W = fliplr(fscanf(fp,'%d %d',[1 2]));
  % read data
  W = fscanf(fp,'%g',size_W);

  if ~feof(fp)
    [size_B,c_B] = fscanf(fp,'%d %d',[1 2]);
    size_B = fliplr(size_B);
    if c_B==2
      assert(~any(size_W));
      assert(isempty(W));
      % Finish reading header: read '\n' char
      [lf,nlf] = fread(fp,1,'char*1');
      assert(nlf==1);
      assert(lf==int8(sprintf('\n')));
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
    error('Size in header (%d,%d) did not match size of data (%d,%d) in file',size_W,size(W));
  end

end
