function B = is_boundary_facet(E,F)
  % IS_BOUNDARY_FACET Determine for each edge E if it is a boundary edge in F
  % for tets. Edges are now facets.
  %
  % B = is_boundary_facet(E,F)
  %
  % Inputs:
  %   E  #E by element-size-1 list of facets
  %   F  #F by element-size list of elements
  % Outputs:
  %   B  #E list bools. true iff unoriented facet occurs exactly once in F
  %     (non-manifold and non-existant edges will be false)
  %

  switch size(F,2)
  case 3
    allE = [F(:,[2 3]);F(:,[3 1]);F(:,[1 2])];
    %% Alternatively (slightly slower)
    %O = outline(F);
    %B = ismember(sort(E,2),sort(O,2),'rows');
  case 4
    allE = [ ...
      F(:,2) F(:,4) F(:,3); ...
      F(:,1) F(:,3) F(:,4); ...
      F(:,1) F(:,4) F(:,2); ...
      F(:,1) F(:,2) F(:,3); ...
      ];
  end
  % http://www.mathworks.com/matlabcentral/newsreader/view_thread/165556
  [~,~,EMAP]=unique(sort([E;allE],2),'rows');
  N = accumarray(EMAP,1);
  % Look of occurances of 2: one for original and another for boundary
  B = N(EMAP(1:size(E,1)))==2;

end
