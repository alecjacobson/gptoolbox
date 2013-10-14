function [I,R,C] = in_elements(F,T)
  % IN_ELEMENTS  Check whether each facet in F truly apears as a facet of the
  % at least one of the elements in T
  %
  % [I] = in_elements(F,T)
  % [I,R,C] = in_elements(F,T)
  %
  % Inputs:
  %   F  #F by dim list of facets
  %   T  #T by dim+1 list of elements
  % Outputs:
  %   I  #F list of indicators whether facet is in element list
  %   R  #F list revealing which *first* tet
  %   C  #F list revealing where in *first* tet
  %

  assert(size(F,2)+1 == size(T,2));

  switch size(F,2)
  case 2
    allF = [T(:,[2 3]);T(:,[3 1]);T(:,[1 2])];
  case 3
    allF = [T(:,[2 3 4]);T(:,[3 4 1]);T(:,[4 1 2]);T(:,[1 2 3])];
  end
  [I,LOCB] = ismember(sort(F,2),sort(allF,2),'rows');
  R = mod(LOCB-1,size(T,1))+1;
  C = floor((LOCB-1)/size(T,1))+1; 

end
