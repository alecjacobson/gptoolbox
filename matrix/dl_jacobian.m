function dfdx = dl_jacobian(obj,x)
  % dfdx = dl_jacobian(obj,x)
  %
  % Inputs:
  %  obj  function handle callable with f = obj(x)
  %  x  #x query point 
  % Outputs:
  %  dfdx  #f by #x autodiff approximation of the jacobian of obj at x
  % 

  dl_x = dlarray(x);
  res = obj(dl_x)
end


