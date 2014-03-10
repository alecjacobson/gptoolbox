function writeTGF(filename,V,E)
  % WRITETGF Write a graph to a .tgf file
  %
  % writeTGF(filename,V,E)
  %
  % Input:
  %  filename  .tgf file name
  %  V  # vertices by 3 list of vertex positions
  %  E  # edges by 2 list of edge indices
  % 
  % Assumes that graph vertices are 3 dimensional

  if( size(V,2) ~= 3)
    % append zeros
    V = [V 0*V(:,1)];
  end
  %disp(['writing: ',filename]);
  fp = fopen(filename,'w');
  % print vertices
  fprintf(fp,'%d %g %g %g\n',[1:size(V,1); V']);
  % print separator
  fprintf(fp,'#\n');
  % print edges, .tgf is 1-indexed
  if(~isempty(E))
    fprintf(fp,'%d %d\n',[E']);
  end
  fclose(fp);
end
