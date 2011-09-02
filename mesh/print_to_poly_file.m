function print_to_poly_file(V, segments, poly_file_name,holes)
  % prints a vertices to a .poly file, with segments connecting those vertices
  %
  % Input
  %   V:               list of vertex positions (2d)
  %   segments:        list of segments (vertex index pairs)
  %   poly_file_name:  name of output file as string (caution! will clobber
  %                    existing)
  %
  error('This should no longer be used. Use writePOLY instead');

  % open file for writing
  poly_file_handle = fopen(poly_file_name,'w');

  % vertices section
  fprintf(poly_file_handle,'# vertices\n');
  fprintf(poly_file_handle,'%d 2 0 0\n', size(V,1));
  for j=1:size(V,1),
    fprintf(poly_file_handle,'%d %.17f %.17f\n',j,V(j,1),V(j,2));
  end

  % segments section
  fprintf(poly_file_handle,'# segments\n');
  fprintf(poly_file_handle,'%d 0\n', size(segments,1));
  for j=1:size(segments,1),
    fprintf(poly_file_handle,'%d %d %d\n',j,segments(j,1),segments(j,2));
  end

  % holes section
  if(exist('holes'))
    fprintf(poly_file_handle,'# holes\n');
    fprintf(poly_file_handle,'%d 0\n', size(holes,1));
    for j=1:size(holes,1),
      fprintf(poly_file_handle,'%d %.17f %.17f\n',j,holes(j,1),holes(j,2));
    end
  else
    % mandatory number of holes lines
    fprintf(poly_file_handle,'# holes\n');
    fprintf(poly_file_handle,'0\n');
  end

  fprintf(poly_file_handle,'\n');
  fclose(poly_file_handle);

end
