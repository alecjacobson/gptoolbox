function s = path_to_triangle()
  % PATH_TO_TRIANGLE Returns absolute, system-dependent path to triangle
  % executable
  %
  % Outputs:
  %   s  path to triangle as string
  %  
  % See also: triangle

  if ispc
    warning([ ...
      'Dear Ladislav, is there a standard place to put executables on a pc?' ...
      'Could you put triangle there and change this accordingly?' ...
      'Thanks, Alec']);
    s = 'c:/prg/lib/triangle/Release/triangle.exe';
  elseif ismac
    s = '/opt/local/bin/triangle';
  elseif isunix 
    % I guess this means linux
    s = '/usr/local/bin/triangle';
  end
end

