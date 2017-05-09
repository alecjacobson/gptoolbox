function str = fphong()
  % FPHONG  Simply returns:
  % struct('FaceColor','interp','FaceLighting','phong');
  %
  % str = fphong()
  % 
  % Example:
  %  set(tsurf(F,V),fphong);
  % 
  str = struct('FaceColor','interp','FaceLighting','phong');
end
