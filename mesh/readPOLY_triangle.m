function [V,F,BM] = readPOLY_triangle(filename)
  % READPOLY_TRIANGLE Read .poly file of triangle output
  %
  % Inputs:
  %   filename  path to .poly file
  % Outputs:
  %   V  #V list of vertices
  %   F  #F by simplex_size list of facets
  %   BM  #F by 1 list of boundary markers
  %   H  #H by dim list of holes (not supported)
  %

  fp = fopen(filename,'r');
  line = eat_comments(fp,'#');
  [Vhead,count] = sscanf(line,'%d %d %d %d');
  assert(count==4);
  [V,count] = fscanf(fp,'%g',[(1+sum(Vhead(2:end))) Vhead(1)]);
  V = V(2:Vhead(2)+1,:)';
  assert(count == (Vhead(1)*(1+sum(Vhead(2:end)))));
  line = eat_comments(fp,'#');
  [Fhead,count] = sscanf(line,'%d %d');
  assert(count == 2);
  [F,count] = fscanf(fp,'%d',[3+Fhead(2) Fhead(1)]);
  % only segments are supported
  BM = F(4:3+Fhead(2),:)';
  F = F(2:3,:)';
  assert(count == prod([3+Fhead(2) Fhead(1)]));
  line = eat_comments(fp,'#');
  [Hhead,count] = sscanf(line,'%d %d');
  assert(Hhead(1) == 0);
  line = eat_comments(fp,'#');
  if ~isempty(line)
    [Rhead,count] = sscanf(line,'%d %d');
    assert(Rhead(1) == 0);
  end

  fclose(fp);

end
