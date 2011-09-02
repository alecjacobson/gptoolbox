function printDiagonal(file_name,D)
  assert(size(D,1)==size(D,2));
  file_id = fopen(file_name,'wt');
  fprintf(file_id,'%d\n', min(size(D)));
  fprintf(file_id,'%.15ld\n',full(diag(D)));
  fclose(file_id);
end
