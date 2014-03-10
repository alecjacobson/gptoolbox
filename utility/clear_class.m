function clear_class(class_name)
  % CLEAR_CLASS  clears all variables of a specific class
  % 
  % clear_class(class_name)
  % 
  % Inputs:
  %   class_name  name of class to clear as string
  %
  % See also: clear
  %

  % Get a cell list of all variable names as strings
  s_v = who;
  s_v
  % loop over variabls
  for ii = 1:numel(s_v)
    % Get this variable as string
    s_s = s_v{ii};
    s_s
    s_d = eval(s_s);
    class(s_d)
    if strcmp(class(s_d),class_name)
      fprintf('clearing: %s...',s_s);
      clear(s_s);
    end
  end
  
end
