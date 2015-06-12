function s = path_to_eltopo()
  % PATH_TO_eltopo Returns absolute, system-dependent path to eltopo header and
  % includes
  %
  % Outputs:
  %   s  path to eltopo base directory as string
  %  
  % See also: eltopo

  if ispc
    warning([ ...
      'Dear Ladislav, is there a standard place to put executables on a pc?' ...
      'Could you put eltopo there and change this accordingly?' ...
      'Thanks, Alec']);
    s = 'c:/prg/lib/eltopo/'
  elseif ismac
    s = '/usr/local/eltopo';
  end
end

