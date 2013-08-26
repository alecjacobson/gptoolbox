function [SV,SVI,SVJ] = remove_duplicate_vertices(V,epsilon)
  % REMOVE_DUPLICATE_VERTICES Remove duplicate vertices upto a uniqueness
  % tolerance (epsilon)
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   epsilon  uniqueness tolerance (significant digit)
  % Outputs:
  %   SV  #SV by dim new list of vertex positions
  %   SVI #V by 1 list of indices so SV = V(SVI,:) 
  %   SVJ #SV by 1 list of indices so V = SV(SVJ,:)
  %
  % Example:
  %   % Mesh in (V,F)
  %   [SV,SVI,SVJ] = remove_duplicate_vertices(V,1e-7);
  %   % remap faces
  %   SF = SVJ(F);
  %
  assert(nargin==1 || epsilon >= 0);
  if nargin==1 || epsilon == 0
    [SV,SVI,SVJ] = unique(V,'rows','stable');
  else
    [~,SVI,SVJ] = unique(round(V/(10*epsilon)),'rows','stable');
    SV = V(SVI,:);
  end
end
