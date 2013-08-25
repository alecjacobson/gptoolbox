function [V,E,F,H] = readPOLY_pyramid(filename)
  % READPOLY_PYRAMID Read .poly file of pyramid/svr output
  %
  % Inputs:
  %   filename  path to .poly file
  % Outputs:
  %   V  #V list of vertices
  %   E  #E by 2+boundary_markers list of segment indices, indexing V, and
  %     optional boundary markers
  %   F  #F struct containing polygon information arrays
  %     .facets  a #facets list of facets,  each facet is a polygon
  %       **NOTE: facets index E *not* V, contrary to typical (V,F) meshes and
  %       contrary to writePOLY_tetgen prototype
  %     .boundary_markers a #facets list of boundary_markers
  %     OR
  %   F  #F by constant-degree+boundary_markers  list of facets
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
  [Ehead,count] = sscanf(line,'%d %d');
  assert(count == 2);
  [E,count] = fscanf(fp,'%d',[3+Ehead(2) Ehead(1)]);
  % only segments are supported
  BME = E(4:3+Ehead(2),:)';
  E = E(2:3,:)';
  offset = 1-min(min(E(:,1)),1);
  if offset
    warning('Offseting indices by %d',offset);
    E = E+offset;
  end
  assert(count == prod([3+Ehead(2) Ehead(1)]));
  line = eat_comments(fp,'#');
  [Fhead,count] = sscanf(line,'%d %d');
  % Not sure if boundary markers are allowed here
  if count == 1
    Fhead(2) = 0;
  end
  % preallocate
  %F = repmat(struct('facets',[],'boundary_markers',[]),Fhead(1),1);
  F.facets = cell(Fhead(1),1);
  F.boundary_markers = cell(Fhead(1),1);
  min_F = inf;
  for f = 1:Fhead(1)
    line = eat_comments(fp,'#');
    [EF,count] = sscanf(line,'%d');
    min_F = min(min_F,EF(1));
    assert(count >= 2);
    assert(count == 2+EF(2));
    F.facets{f} = EF(3:end);
    % not supported
    F.boundary_markers{f} = [];
  end

  fs = cell2mat(cellfun(@size,F.facets,'UniformOutput',false));
  % try to collapse F
  if all(fs(:,1) == 1) && all(fs(:,2) == fs(1,2))
      F = cell2mat(F.facets);
  end

  offset = 1-min(min_F,1);
  if offset
    warning('Offseting indices by %d',offset);
    F.facets = cellfun(@plus,F.facets, ...
      mat2cell(offset*ones(numel(F.facets),1),ones(numel(F.facets),1),1), ...
      'UniformOutput',false);
  end

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
