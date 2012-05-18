function [LF,in] = limit_faces(F,L,exclusive)
  % LIMIT_FACES limit given faces F to those which contain (only) indices found
  % in L.
  %
  % [LF] = limit_faces(F,L,exclusive);
  % [LF,in] = limit_faces(F,L,exclusive);
  %
  % Inputs:
  %   F  #F by 3 list of face indices
  %   L  #L by 1 list of allowed indices
  %    or
  %   L  #V by 1 list of bools specifying if each index in L is allowed
  %   exclusive  flag specifying whether a face is included only if all its
  %     indices are in L, default is false
  % Outputs:
  %   LF  #LF by 3 list of remaining faces after limiting
  %   in  #F list of whether given face was included
  %

  if islogical(L)
    indices = 1:numel(L);
    L = indices(L);
  end


  if(exist('exclusive','var') && exclusive)
    % all indices must be in F
    in = ismember(F(:,1),L) & ismember(F(:,2),L) & ismember(F(:,3),L);
  else
    % at least one index must be in F
    in = ismember(F(:,1),L) | ismember(F(:,2),L) | ismember(F(:,3),L);
  end

  LF = F(in,:);

end
