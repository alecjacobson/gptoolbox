function [SF,SFI,SFJ] = remove_duplicate_simplices(F)
  % REMOVE_DUPLICATE_SIMPLICES Remove duplicate simplices (regardless of order)
  %
  % Inputs:
  %   F  #F by ss list of simplex indices
  % Outputs:
  %   SF  #SF by ss list of unique simplices (maintaining original order of
  %     first simplex found)
  %   SFI  #SF list so that SF = F(SFI,:)
  %   SFJ  #SFI list so that F = SF(SFJ,:)
  [~,SFI,SFJ] = unique(sort(F,2),'rows','stable');
  SF = F(SFI,:);
end
