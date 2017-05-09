function str = fsoft()
  % FSOFT Simply returns:
  % struct('SpecularStrength',0.2,'DiffuseStrength',0.2,'AmbientStrength',0.8);
  %
  % str = fsoft()
  % 
  % Example:
  %  set(tsurf(F,V),fsoft);
  % 
  str = struct('SpecularStrength',0.2,'DiffuseStrength',0.2,'AmbientStrength',0.8);
end

