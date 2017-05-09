function writePOLY_pyramid(filename,V,E,F,H)
  % WRITEPOLY_PYRAMID prints vertices, segments and facets to a .poly suitable
  % for use with PYRAMID
  %
  % writePOLY_pyramid(filename,V,E,F,H)
  %
  % Inputs:
  %   V  #V by dim=3 list of vertex positions
  %   E  #E by 2+boundary_markers list of segment indices, indexing V, and
  %     optional boundary markers
  %   F  #F struct containing polygon information arrays
  %     .facets  a #facets list of facets,  each facet is a polygon
  %       **NOTE: facets index E *not* V, contrary to typical (V,F) meshes and
  %       contrary to writePOLY_tetgen prototype
  %     .boundary_markers a #facets list of boundary_markers
  %     OR
  %   F  #F by constant-degree+boundary_markers  list of facets
  %
  % Example:
  %   % Mesh in (V,F)
  %   % Gather all edges
  %   E = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
  %   [~,IE,IuE] = unique(sort(E,2),'rows');
  %   % unique edges
  %   uE = E(IE,:);
  %   % reindex F into E
  %   FinE = reshape(IuE,size(F));
  %   Facets = [];
  %   Facets.facets = mat2cell(FinE,ones(size(FinE,1),1),[3]);
  %   Facets.boundary_marker = -ones(size(FinE,1),1);
  %   writePOLY_pyramid(path,V,uE,Facets,[]);
  %
  % See also: cdt, tetgen, writePOLY_triangle, writePOLY_tetgen
  %
  %

  % open file for writing
  poly_file_handle = fopen(filename,'w');

  dim = size(V,2);


  if dim ~= 3
    error('writePOLY_pyramid is for 3d meshes. Try writePOLY_triangle etc.');
  end

  % vertices section
  fprintf(poly_file_handle,'# Part 1 - node list\n');
  format = '%d %.17g %.17g %.17g\n';
  fprintf(poly_file_handle,'%d %d 0 0\n', size(V,1),size(V,2));
  if ~isempty(V)
    fprintf(poly_file_handle,format,[1:size(V,1);V']);
  end

  fprintf(poly_file_handle,'# Part 2 - segment list\n');
  fprintf(poly_file_handle,'%d %d\n',size(E,1),size(E,2)-2);
  format = ['%d %d %d' repmat(' %d',1,size(E,2)-2) '\n'];
  fprintf(poly_file_handle,format, [1:size(E,1);E']);

  fprintf(poly_file_handle,'# Part 2 - facet list\n');
  % for now, always include boundary markers
  % [num facets] [boundary markers]
  assert(isempty(F) || isempty(F.facets) || iscell(F.facets));
  if isempty(F)
    fprintf(poly_file_handle,'0\n');
  else
      fprintf(poly_file_handle,'%d %d\n',numel(F.facets),~isempty(F.boundary_marker));
      fs = cell2mat(cellfun(@size,F.facets,'UniformOutput',false));
      % Try to print all at once if facets are all the same size
      if ~isempty(fs) && all(fs(:,1) == 1) && all(fs(:,2) == fs(1,2))
          % build format
          fformat = [ ...
              ... % index and size
              '%d ' num2str(fs(1,2)) ...
              ... % facet indices into segments
              repmat(' %d',1,fs(1,2)) ...
              ... % boundary markers
              repmat(' %d',1,size(F.boundary_marker,2)) '\n'];
          % print all at once
          fprintf(poly_file_handle,fformat, ...
              [1:size(F.facets,1);[cell2mat(F.facets) F.boundary_marker]']);
      else
          % irregular face valences
          for f=1:numel(F.facets)
              % 1d list
              assert(size(F.facets{f},1)==1 || size(F.facets{f},2) == 1);
              fprintf('%d %d',f, numel(F.facets{f}));
              % print indices
              for p=1:numel(F.facets{f})
                  fprintf(poly_file_handle,' %d',F.facets{f}(p));
              end
              fprintf(poly_file_handle,'\n');
          end
      end
  end

  % [num holes]
  fprintf(poly_file_handle,'# Part 3 - hole list\n');
  fprintf(poly_file_handle,'%d\n',size(H,1));
  if ~isempty(H)
    assert(isempty(V) || size(H,2) == size(V,2));
    fprintf(poly_file_handle,'%0.17g %0.17g %0.17g\n',[1:size(H,1);H']);
  end

end
