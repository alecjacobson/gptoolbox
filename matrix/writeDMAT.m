function writeDMAT(filename,W,ascii)
  % WRITEDMAT  writes a matrix from a dmat file.  first line is <# columns> <#
  % rows>, then values with columns running faster
  %
  % writeDMAT(filename,W)
  % writeDMAT(filename,W,asci)
  %
  % Input:
  %   filename  name of .dmat file
  %   W  matrix read from file
  %   ascii  write ascii file {true}
  %
  % Example:
  %   % weights for bones (C,BE), size #V by #BE
  %   WB;
  %   % rearrange weights to match sorted bone edge order used by lbs.app
  %   [~,I] = sortrows(BE,[1 2]);
  %   WB = WB(:,I);
  %   % Append zeros for dummy point handle weights
  %   W = [zeros(size(V,1),size(C,1)) WB];
  %   % finally, write weights to file, usable by lbs.app
  %   writeDMAT(filename,W);
  %
  % See also: readDMAT
  %
  
  %disp(['writing: ',filename]);
    % open file
    fp = fopen(filename,'w');
  if ~exist('ascii','var') || ascii
    % write header
    fprintf(fp,'%d %d\n',fliplr(size(W)));
    % write data
    fprintf(fp,'%.16g\n',W);
  else
    % write header for ascii
    fprintf(fp,'0 0\n');
    % write header for binary
    fprintf(fp,'%d %d\n',fliplr(size(W)));
    if ~isa(W,'double')
      warning('Input will be cast from %s to double',class(W));
    end
    fwrite(fp,double(W),'double');
  end
  % close file
  fclose(fp);

end
