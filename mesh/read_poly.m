function poly = read_poly(poly_file_name)
  % reads vertices to a .poly file, with segments connecting those vertices
  %
  % Input
  %   poly_file_name:  name of output file as string (caution! will clobber
  %                    existing)
  % Output
  %   poly:      struct array where each element contains fields:
  %                    x,y,hole
  %

  % open file for reading
  poly_file_handle = fopen(poly_file_name,'r');
  poly = [];

  % read vertices
  % read segments
  % read holes

  fclose(poly_file_handle);
end
