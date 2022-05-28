function E = material_boundary_facets(F,FI)
  % E = material_boundary_facets(F,FI)
  %
  % Inputs:
  %   F  #F by 3 list of triangle indices into rows of some V
  %   FI #F list of material ids
  % Outputs:
  %   E  #E by 2 list of unique edge boundaries between elements with different
  %     ids
  [uE,~,EMAP] = unique(sort([F(:,[2 3]);F(:,[3 1]);F(:,[1 2])],2),'rows');
  uE2FI = sparse(EMAP,repmat(FI,3,1),1);
  E = uE(sum(uE2FI,2)==1 | sum(uE2FI>0,2)==2,:);
end
