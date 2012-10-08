function writeBB(filename,D,type)
  % WRITEBB  write medit data file
  %
  % writeBB(filename,D,type)
  %
  % Inputs:
  %   filename  path to *.bb file
  %   D  #D by dim  data matrix
  %   type  type identifier as per
  %     http://www.ann.jussieu.fr/frey/publications/RT-0253.pdf page 34
  %     2 means defined at vertices
  %     1 means defined at elements

  f = fopen(filename,'w');
  % print header
  % only vector data is supported
  fprintf(f,'%d %d %d %d\n',2,1,size(D,1),type);
  fprintf(f,'%g\n',D');
  fclose(f);

end
