function writeELE(filename,E)
  % WRITEELE Write mesh elements to a .ele file
  %
  % writeELE(filename,V)
  %
  % Inputs:
  %  filename  name of output file
  %  E  list of elements
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %

  fp = fopen(filename,'w');
  % attributes are not supported
  % number of edges number of elements  0 0
  fprintf(fp,'%d %d 0\n',size(E));
  % .node is 1-indexed
  indices = 1:size(E,1);
  % build format string
  str = '%d';
  for(ii = 1:size(E,2))
    str = [str ' %d'];
  end
  str = [str '\n'];
  fprintf(fp,str,[indices', E]');
  fclose(fp);
end
