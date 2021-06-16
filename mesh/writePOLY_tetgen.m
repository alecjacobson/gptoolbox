function writePOLY_tetgen(filename,V,F,H,varargin)
  % WRITEPOLY_TETGEN prints vertices and planar facets to a .poly file for
  % tetgen
  %
  % writePOLY_tetgen(filename,V,F,H)
  %
  % Input
  %   filename:  name of output file as string (caution! will clobber
  %                    existing)
  %   V  #V by dim list of vertex positions
  %   F  #F struct containing polygon information arrays
  %     .facets  a #facets list of facets, each facet is a again a list of
  %       polygons
  %      ** Note: contrary to writePOLY_pyramid, here facets index V directly
  %     .boundary_markers a #facets list of boundary_markers
  %     .holes  a #facets list of holes, each holes entry is again a list for
  %       each facet
  %    or 
  %   F  #F by uniform_facet_size list of facets into V
  %   H  #H by dim list of volume hole positions
  %   Optional:
  %     'BoundaryMarkers' followed by #facets list of boundary markers
  %     'BoundaryMarkersNodes' followed by #facets list of boundary markers
  %     'RegionList' followed by #regions list of region vertex positions and (optional) region numbers and attributes
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

  v = 1;
  BM = [];
  BMN = [];
  RL = [];
  while v <= numel(varargin)
    switch varargin{v}
    case 'BoundaryMarkers'
      assert((v+1)<=numel(varargin));
      v = v+1;
      BM = varargin{v};
    case 'BoundaryMarkersNodes'
      assert((v+1)<=numel(varargin));
      v = v+1;
      BMN = varargin{v};
    case 'RegionList'
      assert((v+1)<=numel(varargin));
      v = v+1;
      RL = varargin{v};
    otherwise
      error(['Unsupported parameter: ' varargin{v}]);
    end
    v=v+1;
  end

  % open file for writing
  poly_file_handle = fopen(filename,'w');

  % dimensions in V, should be 3
  dim = size(V,2);
  if isempty(V)
    dim = 3;
  end

  if dim ~= 3
    error('writePOLY_tetgen is for 3d meshes. Try writePOLY_triangle etc.');
  end

  % vertices section
  fprintf(poly_file_handle,'# vertices\n');
  fprintf(poly_file_handle,'# Part 1 - node list\n');
  fprintf(poly_file_handle,'%d %d 0 %d\n', size(V,1),size(V,2),~isempty(BMN));
  format = '%d %.17g %.17g %.17g\n';
  if ~isempty(V)
      if ~isempty(BMN)
          assert(numel(BMN)==size(V,1));
          formatV = '%d %.17g %.17g %.17g %d\n';
          fprintf(poly_file_handle,formatV,[1:size(V,1);V';BMN']);
      else
          formatV=format;
          fprintf(poly_file_handle,formatV,[1:size(V,1);V']);
          
      end
  end

  fprintf(poly_file_handle,'# Part 2 - facet list\n');
  % Try to print all at once if facets are all the same size
  if ~isstruct(F)
    fprintf(poly_file_handle,'%d %d\n',size(F,1),~isempty(BM));
    assert(numel(BM)==size(F,1) || isempty(BM));
    % build format
    fformat = '1 0';
    if ~isempty(BM)
      fformat = [fformat ' %d'];
    end
    fformat = [fformat '\n' num2str(size(F,2))];
    for p=1:size(F,2)
      fformat = [fformat ' %d']; %#ok<AGROW>
    end
    fformat = [fformat '\n'];
    % print all at once
    fprintf(poly_file_handle,fformat,[BM F]');
  else
    fprintf(poly_file_handle,'%d %d\n',numel(F.facets),1);
    % irregular face valences
    for f=1:numel(F.facets)
      % [num polygons] [num holes] [boundary marker]
      fprintf(poly_file_handle,'%d %d %d\n', ...
        numel(F.facets{f}),size(F.holes{f},1),F.boundary_marker(f));
      % loop over polygons
      for p=1:numel(F.facets{f})
        % [num corners] [corner 1] [corner 2] ...
        fprintf(poly_file_handle,' %d',numel(F.facets{f}{p}));
        fprintf(poly_file_handle,' %d',F.facets{f}{p});
        fprintf(poly_file_handle,'\n');
      end
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
  % [num regions]
  fprintf(poly_file_handle,'# Part 4 - region list\n');
  fprintf(poly_file_handle,'%d\n',size(RL,1));
  if ~isempty(RL)
    % [region #] [region x] [region y] [region z] [region attribute] and/or [region number]
    formatRL = ['%d' repmat(' %.17g',1,size(RL,2)) '\n'];
    fprintf(poly_file_handle,formatRL,[1:size(RL,1);RL']);
  end
  fprintf(poly_file_handle,'\n');
  fclose(poly_file_handle);
end
