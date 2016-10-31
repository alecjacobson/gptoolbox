function [V,E,BME] = readPOLY_triangle(filename)
  % READPOLY_TRIANGLE Read .poly file of triangle output
  %
  % Inputs:
  %   filename  path to .poly file
  % Outputs:
  %   V  #V list of vertices
  %   E  #E by 2 list of segments
  %   BME  #E by 1 list of boundary markers
  %   H  #H by dim list of holes (not supported)
  %
  % Example:
  %   % ... triangle on (V,E)
  %   [VV,EE,BME] = readPOLY_triangle(...)
  %   % Triangle does not maintain orders of edges 
  %   EEdotE = sum((VV(EE(:,1),:)-VV(EE(:,2),:)).*(V(E(BME,1),:)-V(E(BME,2),:)),2);
  %   EE = bsxfun(@times,(EEdotE>0),EE) + bsxfun(@times,(EEdotE<=0),fliplr(EE));

  fp = fopen(filename,'r');
  line = eat_comments(fp,'#');
  [Vhead,count] = sscanf(line,'%d %d %d %d');
  assert(count==4);
  [V,count] = fscanf(fp,'%g',[(1+sum(Vhead(2:end))) Vhead(1)]);
  V = V(2:Vhead(2)+1,:)';
  assert(count == (Vhead(1)*(1+sum(Vhead(2:end)))));
  line = eat_comments(fp,'#');
  [Ehead,count] = sscanf(line,'%d %d');
  assert(count == 2);
  [E,count] = fscanf(fp,'%d',[3+Ehead(2) Ehead(1)]);
  % only segments are supported
  BME = E(4:3+Ehead(2),:)';
  E = E(2:3,:)';
  assert(count == prod([3+Ehead(2) Ehead(1)]));
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
