function [S,uF,I,J] = total_signed_occurrences(F)
  % TOTAL_SIGNED_OCCURRENCES Compute the total number of times each polygon in F
  % occurs.
  %
  % Inputs:
  %   F  #F by d list of simple polygon indices
  % Outputs:
  %   S  #uF list of total signed occurances
  %   uF  #uF by d list of unique polygons (sorted)
  %   I  #uF list of indices into F
  %   J  #F list of indices into uF
  %
  m = size(F,1);
  d = size(F,2);
  % In general d>3 can't just sort
  sF = F;
  MN = min(sF,[],2);
  % rotate until smallest index is first
  for p = 1:d
    bad = MN~=sF(:,1);
    sF(bad,:) = sF(bad,[2:d 1]);
  end
  bad = sF(:,d)<sF(:,2);
  sF(bad,2:d) = fliplr(sF(bad,2:d));

  % matches sort modulo rotation?
  M = zeros(m,1);
  P = [1:d];
  for p = 1:d
    M = M | all(sF(:,P) == F,2);
    % rotate
    P = [P(2:end) P(1)];
  end


  [uF,I,J] = unique(sF,'rows','stable');
  S = full(sparse(J,1,2*M-1,size(uF,1),1));
end
