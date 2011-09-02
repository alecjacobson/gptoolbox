function printIJV(file_name,S)
  file_id = fopen(file_name,'wt');
  fprintf(file_id,'%d %d\n', size(S));
  [i,j,v] = find(S);
  fprintf(file_id,'%d %d %.15ld\n',[i-1,j-1,v]');
  fclose(file_id);
end
