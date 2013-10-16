function F = triangles_from_strip(S)
  % TRIANGLES_FROM_STRIP
  %
  % Inputs:
  %   S  #S list of indices
  % Outputs:
  %   F  #S-2 by 3 list of triangle indices
  %

  F = zeros(numel(S)-2,3);
  for s = 3:numel(S)
    F(s-2,:) = [S(s-2) S(s-1) S(s)];
  end

end
