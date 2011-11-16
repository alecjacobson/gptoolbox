function writeDMAT(filename,W)
  % WRITEDMAT  writes a matrix from a dmat file.
  %   first line is <# columns> <# rows>, then values with columns running
  %   faster
  %
  % writeDMAT(filename,W)
  %
  % Input:
  %   filename  name of .dmat file
  %   W  matrix read from file
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
  
  disp(['writing: ',filename]);
  % open file
  fp = fopen(filename,'w');
  % write header
  fprintf(fp,'%d %d\n',fliplr(size(W)));
  % write data
  fprintf(fp,'%.16g\n',W);
  % close file
  fclose(fp);

end
