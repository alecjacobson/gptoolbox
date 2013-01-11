function writePOLY_tetgen(varargin)
  % WRITEPOLY prints a vertices to a .poly file, with E connecting those
  % vertices or planar facets F
  %
  % writePOLY_tetgen(poly_file_name,poly_struct)
  % writePOLY_tetgen(poly_file_name,V,E,H)
  % writePOLY_tetgen(poly_file_name,V,F,H)
  %
  % Input
  %   poly_file_name:  name of output file as string (caution! will clobber
  %                    existing)
  %   poly:      struct array where each element contains fields:
  %                    x,y,hole
  %     OR
  %   V  #V by dim=3 list of vertex positions
  %   E  #E by 2 list of edges E
  %   H  #H by dim list of hole positions
  %     OR
  %   V  #V by dim list of vertex positions
  %   F  #F struct containing polygon information arrays
  %     .facets  a #facets list of facets, each facet is a again a list of
  %       polygons (currently all polygons of a facet must be same valence)
  %      ** Note: contrary to writePOLY_pyramid, here facets index V directly
  %     .boundary_markers a #facets list of boundary_markers
  %     .holes  a #facets list of holes, each holes entry is again a list for
  %       each facet
  %   H  #H by dim list of hole positions
  %
  % Example:
  %   % constrained Delaunay tetrahedralization of a triangle mesh
  %   [V,F] = load_mesh('gargoyle.obj');
  %   % constrained Delaunay tetrahedralization of just points
  %   D = DelaunayTri(V);
  %   % Faces on boundary of convex hull
  %   BF = boundary_faces(D.Triangulation);
  %   % avoid dupicate faces
  %   FmBF = setdiff(sort(F,2),sort(BF,2),'rows');
  %   Facets = [];
  %   Facets.facets = mat2cell([BF;FmBF],ones(size(BF,1)+size(FmBF,1),1),[3]);
  %   Facets.boundary_marker = [ones(size(BF,1),1);-ones(size(FmBF,1),1)];
  %   Facets.holes = cell(numel(Facets.facets),1);
  %   writePOLY_tetgen('temp.poly',V,Facets,[]);
  %   % -p: we're giving PLC, -g: output .mesh, -Y: no steiners
  %   !/usr/local/bin/tetgen -pgY ~/Documents/volume/temp.poly
  %   [VV,TT,FF] = readMESH('temp.1.mesh');
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: cdt, tetgen, writePOLY_triangle, writePOLY_pyramid
  %

  if(nargin == 4)
    poly_file_name = varargin{1};
    V = varargin{2};
    E = varargin{3};
    H = varargin{4};
  elseif(nargin == 2)
    poly_file_name = varargin{1};
    poly = varargin{2};
    [V,E,H] = poly2VEH(poly);
  else
    error('Wrong number of inputs');
  end
  % open file for writing
  poly_file_handle = fopen(poly_file_name,'w');

  % dimensions in V, should be 3
  dim = size(V,2);
  if isempty(V)
    if isstruct(E)
      % min polygon size
      ms = min(cell2mat(cellfun(@size,E.facets,'UniformOutput',false)),[],1);
      dim = ms(2);
    else
      dim = size(E,2);
    end
  end

  if dim ~= 3
    error('writePOLY_tetgen is for 3d meshes. Try writePOLY_triangle etc.');
  end

  % vertices section
  fprintf(poly_file_handle,'# vertices\n');
  format = '%d %.17g %.17g %.17g\n';

  fprintf(poly_file_handle,'# Part 1 - node list\n');
  fprintf(poly_file_handle,'%d %d 0 0\n', size(V,1),size(V,2));
  if ~isempty(V)
    fprintf(poly_file_handle,format,[1:size(V,1);V']);
  end
  if isstruct(E)
    F = E;
  else
    F = [];
    F.facets = mat2cell(E,ones(size(E,1),1),size(E,2));
    F.boundary_marker = ones(size(E,1),1);
    F.holes = cell(numel(F.facets),1);
  end
  fprintf(poly_file_handle,'# Part 2 - facet list\n');
  % for now, always include boundary markers
  % [num facets] [boundary markers]
  assert(isempty(F.facets) || iscell(F.facets));
  fprintf(poly_file_handle,'%d %d\n',numel(F.facets),1);
  fs = cell2mat(cellfun(@size,F.facets,'UniformOutput',false));
  fhs = cell2mat(cellfun(@size,F.holes,'UniformOutput',false));
  % Try to print all at once if facets are all the same size
  if ~isempty(fs) && all(fs(:,1) == 1) && all(fs(:,2) == fs(1,2)) && all(fhs(:,1) == 0)
    % build format
    fformat = ['1 0 %d\n ' num2str(fs(1,2))];
    for p=1:fs(1,2)
      fformat = [fformat ' %d'];
    end
    fformat = [fformat '\n'];
    % print all at once
    fprintf(poly_file_handle,fformat, ...
      [F.boundary_marker cell2mat(F.facets)]');
  else
    % irregular face valences
    for f=1:numel(F.facets)
      % [num polygons] [num holes] [boundary marker]
      fprintf(poly_file_handle,'%d %d %d\n', ...
        size(F.facets{f},1),size(F.holes{f},1),F.boundary_marker(f));
      % [num corners] [corner 1] [corner 2] ...
      fprintf(poly_file_handle,'%d',size(F.facets{f},2));
      for p=1:numel(F.facets{f})
        fprintf(poly_file_handle,' %d',F.facets{f}(p));
      end
      fprintf(poly_file_handle,'\n');
      % [hole #] [hole x] [hole y] [hole z]
      if ~isempty(F.holes{f})
        assert(size(F.holes{f},2) == size(V,2));
        fprintf(poly_file_handle,format,[1:size(F.holes{f},1);F.holes{f}']);
      end
    end
  end

  % [num holes]
  fprintf(poly_file_handle,'# Part 3 - hole list\n');
  fprintf(poly_file_handle,'%d\n',size(H,1));
  if ~isempty(H)
    assert(isempty(V) || size(H,2) == size(V,2));
    fprintf(poly_file_handle,format,[1:size(H,1);H']);
  end
  % not supported
  fprintf(poly_file_handle,'# Part 4 - region list\n');
  fprintf(poly_file_handle,'0\n');

  fprintf(poly_file_handle,'\n');
  fclose(poly_file_handle);
end

