function O = bone_offsets(C,P)
  % BONE_OFFSETS  given a list of joint positions and parent indices compute
  % for each joint an offset based on the parent's position (roots will be
  % offset by origin)
  %
  % O = bone_offsets(C,P)
  % O = bone_offsets(C,BE)
  %
  % Inputs:
  %   C  #C by dim list of joint positions
  %   P  #C list of parent indices into C (-1 specifies a root joint)
  %     or 
  %   BE  #B by 2 list of bone indices into C (assumes BE(:,1) contains parents
  %     and BE(:,2) contains children)
  % Outputs:
  %   O  #C by dim list of joint offsets
  %
  % See also: readTGF, writeBF
  %

  % number of joints
  m = size(C,1);

  % second argument is bone list
  if size(P,2) == 2
    if size(P,1) == 1
      warning('Ambiguous varargin{2}');
    end
    BE = P;
    % convert to parent list
    P = -zeros(size(C,1),1);
    P(BE(:,2)) = BE(:,1);
  end

  assert(prod(size(P)) == max(size(P)));
  assert(min(P) == 0);
  P = P(:);
  assert(m == numel(P));

  O = C;
  indices = 1:m;
  % find all non roots
  non_roots = indices(P~=0);
  O(non_roots,:) = C(non_roots,:) - C(P(non_roots),:);

end
