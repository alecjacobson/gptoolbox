function [S,V,res,h] = readSDF(filename)
  % READSDF Read a signed distance field from a [.sdf
  % file](https://github.com/christopherbatty/SDFGen/blob/master/main.cpp#L26) 
  %
  %[S,V,res,h] = readSDF(filename)
  %
  % Input:
  %   filename  path to .sdf file
  % Outputs:
  %   S  #V list of signed distance values 
  %   V  #V by 3 list of grid center locations
  %   res  3-long integer size of grid
  %   h  grid spacing
  %
  fp = fopen( filename, 'r' );
  res = fscanf(fp,'%d %d %d\n',3)';
  origin = fscanf(fp,'%g %g %g\n',3)';
  h = fscanf(fp,'%g\n',1);
  S = fscanf(fp,'%g\n',prod(res));
  assert(size(S,1) == prod(res));
  [X,Y,Z] = meshgrid( ...
    h*linspace(0,res(1)-1,res(1)), ...
    h*linspace(0,res(2)-1,res(2)), ...
    h*linspace(0,res(3)-1,res(3)));
  V = bsxfun(@plus,[X(:) Y(:) Z(:)],origin);
  S = reshape(permute(reshape(S,res),[2 1 3]),[],1);
end
