function s = path_to_eltopo()
  % PATH_TO_eltopo Returns absolute, system-dependent path to eltopo header and
  % includes
  %
  % Outputs:
  %   s  path to eltopo base directory as string
  %  
  % See also: eltopo

  if ispc
    s = 'c:/prg/lib/eltopo/'
  elseif ismac
    s = find_first_path({'/usr/local/eltopo'});
  end

end

