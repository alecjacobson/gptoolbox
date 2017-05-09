function writePOLY_triangle(varargin)
  % WRITEPOLY_TRIANGLE prints a vertices to a .poly file, with E connecting
  % those vertices OR triangels made from F. Suitable for use with TRIANGLE
  %
  % writePOLY_triangle(filename,poly_struct)
  % writePOLY_triangle(filename,V,E,H)
  %
  % Input
  %   filename:  name of output file as string (caution! will clobber
  %                    existing)
  %   poly:      struct array where each element contains fields:
  %                    x,y,hole
  %     OR
  %   V  #V by dim=2 list of vertex positions
  %   E  #E by 2 list of edges E
  %   H  #H by dim list of hole positions
  %     OR
  %   V  #V by dim list of vertex positions
  %   F  #F struct containing polygon information arrays
  %     .facets  a #facets list of facets, each facet is a again a list of
  %       polygons (currently all polygons of a facet must be same valence)
  %     .boundary_markers a #facets list of boundary_markers
  %     .holes  a #facets list of holes, each holes entry is again a list for
  %       each facet
  %   H  #H by dim list of hole positions
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: png2objandtga, png2poly, writePOLY_tetgen
  %

  if(nargin == 4)
    filename = varargin{1};
    V = varargin{2};
    E = varargin{3};
    H = varargin{4};
  elseif(nargin == 2)
    filename = varargin{1};
    poly = varargin{2};
    [V,E,H] = poly2VEH(poly);
  else
    error('Wrong number of inputs');
  end
  % open file for writing
  poly_file_handle = fopen(filename,'w');

  % dimensions in V, should be 2: triangle format
  dim = size(V,2);
  if isempty(V)
    if isstruct(E)
      % min polygon size
      ms = min(cell2mat(cellfun(@size,E.facets,'UniformOutput',false)),[],1);
      dim = ms(2);
    else
      dim = size(E,2);
    end
    if isempty(E)
        dim = 2;
    end
  end

  if dim ~= 2
    error('writePOLY_triangle is for 2d meshes. Try writePOLY_tetgen etc.');
  end

  % vertices section
  format = '%d %.17g %.17g\n';

  fprintf(poly_file_handle,'# Part 1 - node list\n');
  fprintf(poly_file_handle,'%d %d 0 0\n', size(V,1),size(V,2));
  if ~isempty(V)
    fprintf(poly_file_handle,format,[1:size(V,1);V']);
  end
  % The following became convoluted because of ancient merge and unmerge with
  % writePOLY_tetgen.m
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
  % the triangle .poly format is slightly different
  fprintf(poly_file_handle,'%d %d %d %d\n', ...
    [1:numel(F.facets); [cell2mat(F.facets) F.boundary_marker]']);
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

