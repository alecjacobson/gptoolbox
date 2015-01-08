function create_cache(cache_name)
  % CREATE_CACHE Save the workspace of the caller function using the provided
  % cache_name, possibly overwriting it.
  %
  % create_cache(cache_name)
  %
  % Inputs:
  %   cache_name  path to .mat cache file
  %
  % See also: find_cache

  if isunix
    % get a list of current variables in this scope, this is the input "state"
    variables = evalin('caller','who');
    variable_names = sprintf('^%s$|',variables{:});
    save_cmd = sprintf('save(''%s'',''-regexp'',''%s'');',cache_name,variable_names);
    evalin('caller',save_cmd);
  end
end
