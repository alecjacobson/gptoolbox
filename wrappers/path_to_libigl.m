function s = path_to_libigl()
  % PATH_TO_LIBIGL
  %
  % s = path_to_libigl()
  %
  % Returns path to libigl as string
  
  s = find_first_path({'/usr/local/igl/libigl/','/usr/local/libigl/'});
end
