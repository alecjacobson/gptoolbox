function [I,C] = on_boundary(F)
  % ON_BOUNDARY  determines if each face of a manifold mesh is on the boundary
  % (contains at least on boundary edge)
  %
  % [I,C] = on_boundary(F)
  %
  % Inputs:
  %   F  #F by simplex-size list of element indices
  % Outputs:
  %   I  #F list of bools, whether on boundary
  %   C  #F by simplex size matrix of bools, whether opposite edge on boundary
  simplex_size = size(F,2);
  switch simplex_size
  case 3
    E = [F(:,2) F(:,3); F(:,3) F(:,1); F(:,1) F(:,2)];
    [sortedE] = sort(E,2);
    if min(sortedE(:))>0 
      [~,~,n] = unique(sortedE(:,1)+(max(sortedE(:,1)))*sortedE(:,2));
    else
      [~,~,n] = unique(sortedE,'rows');
    end
    counts = accumarray(n(:), 1);
    C = counts(n);
    C = reshape(C,size(F));
    C = C==1;
    I = any(C,2);
  case 4
    T = F;
    allF = [ ...
      T(:,2) T(:,4) T(:,3); ...
      T(:,1) T(:,3) T(:,4); ...
      T(:,1) T(:,4) T(:,2); ...
      T(:,1) T(:,2) T(:,3); ...
      ];
    % sort rows so that faces are reorder in ascending order of indices
    sortedF = sort(allF,2);
    % determine uniqueness of faces
    [u,m,n] = unique(sortedF,'rows');
    % determine counts for each unique face
    counts = accumarray(n(:), 1);
    C = counts(n);
    C = reshape(C,size(T));
    C = C==1;
    I = any(C,2);
  otherwise
    error(['Unsupported simplex size' num2str(simplex_size)]);
  end
end
