function flag = is_writable(name)
  % IS_WRITABLE Determine if a path is writable
  %
  % flag = is_writable(name)
  %
  % Inputs:
  %   name  path to check
  % Outputs:
  %   flag  whether writable
  if exist(name,'file')
    [~,attribs] = fileattrib(name);
    flag = attribs.UserWrite;
    return;
  end
  dirname = fileparts(name);
  if exist(dirname,'file')
    [~,attribs] = fileattrib(dirname);
    flag = attribs.UserWrite;
    return;
  end
  flag = false;
end
