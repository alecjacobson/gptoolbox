function O = occurrence(A,varargin)
  % Mark occurrences of elements in a vector A
  %
  % O = occurrence(A)
  % O = occurrence(A,'rows')
  % 
  % Inputs:
  %   A  #A list of elements
  %   Optional:
  %      'rows'  find occurencess of rows
  % Outputs:
  %   O  #A list of elements so that O(i)=j implies that A(i) is the jth time
  %   we've seen A(i) in A
  if nargin>1 && strcmp(varargin{1},'rows')
    [S,K] = sortrows(A);
    [~,J,I] = unique(S,'stable','rows');
  else
    if ~isvector(A)
      A = A(:);
    end
    [S,K] = sort(A);
    [~,J,I] = unique(S,'stable');
  end
  O = ((1:numel(K))-J(I)')+1;
  O(K) = O;
  if nargin>1 && strcmp(varargin{1},'rows')
    O = reshape(O,size(A,1),1);
  else
    O = reshape(O,size(A));
  end

end
