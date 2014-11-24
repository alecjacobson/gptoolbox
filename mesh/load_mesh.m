function [V,F] = load_mesh(filename,varargin)
  % read in vertices and faces from a .off or .obj file
  % 
  % [V,F] = load_mesh(filename)
  % [V,F] = load_mesh(filename,'ParameterName',ParameterValue, ...)
  %
  % Input:
  %   filename  file holding mesh
  %   Optional:
  %     'Quiet' followed by whether to be quiet {false}
  % Output:
  %   V (vertex list) 
  %   F (face list) fields
  %
  % Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
  %
  % See also: readOBJ, readOBJfast, readOFF
  %

  % parse any optional input
  v = 1;
  quiet = false;
  while(v <= numel(varargin))
    switch varargin{v}
    case 'Quiet'
      v = v+1;
      assert(v<=nargin);
      quiet = varargin{v};
      quiet = quiet * 1;
    otherwise
      error(['Unsupported parameter: ' varargin{v}]);
    end
    v = v + 1;
  end

  [~,~,ext] = fileparts(filename);
  ext = lower(ext);
  switch ext
  case '.off'
    [V,F] = readOFF(filename);
  case '.ply'
    [V,F] = readPLY(filename);
  case '.stl'
    [V,F] = readSTL(filename);
  case '.wrl'
    [V,F] = readWRL(filename);
  case '.obj'
    try
      [V,F] = readOBJfast(filename);
    catch exception
      if ~quiet
        fprintf('Fast reader failed, retrying with more robust, slower reader\n');
      end
      [V,F] = readOBJ(filename);
    end
  otherwise
    error('Unknown mesh format: %s',ext);
  end
end
