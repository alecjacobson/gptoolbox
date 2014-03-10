function F = triangles_from_strip(S)
  % TRIANGLES_FROM_STRIP Create a list of triangles from a stream of indices
  % along a strip.
  %
  % Inputs:
  %   S  #S list of indices
  % Outputs:
  %   F  #S-2 by 3 list of triangle indices
  %

  %F = zeros(numel(S)-2,3);
  %for s = 3:numel(S)
  %  if mod(s,2) == 0
  %    F(s-2,:) = fliplr([S(s-2) S(s-1) S(s)]);
  %  else
  %    F(s-2,:) = [S(s-2) S(s-1) S(s)];
  %  end
  %end

  S = S(:);
  % Maintain order and orientations
  F = [ ...
    S(1:2:end-2) S(2:2:end-1) S(3:2:end) ...
    S(3:2:end-1) S(2:2:end-2) S(4:2:end)];
  F = reshape(F',3,numel(S)-2)';

end
