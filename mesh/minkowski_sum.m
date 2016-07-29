function [SV,SF,J] = minkowski_sum(VA,FA,VB,FB)
  % MINKOWSKI_SUM A very slow and non robust implementation of the Minkowski
  % sum of two triangle meshes. This will likely result in self-intersections
  % and poor performance for even absurdly simple inputs.
  %
  % Inputs:
  %   VA  #VA by 3 list of mesh vertices
  %   FA  #FA by 3 list of triangle indices into VA
  %   VB  #VB by 3 list of mesh vertices
  %   FB  #FB by 3 list of triangle indices into VB
  % Outputs:
  %   SV  #SV by 3 list of mesh vertices
  %   SF  #SF by 3 list of triangle indices into SV
  %   J  #SF list of birth triangles into [FA;edge of B;FB;edge of A]
  % 

  % This is crazy O(n²) slow.

  function [SV,SF,J] = linear_sweep_edges(VA,FA,VB,EB)
    assert(size(EB,2) == 2,'EB must be list of edge indices');
    SV = [];
    SF = [];
    J = [];
    for e = 1:size(EB,1)
      se = VB(EB(e,1),:);
      de = VB(EB(e,2),:);
      v = de-se;
      [SVe,SFe,Je] = linear_sweep(bsxfun(@plus,VA,se),FA,v,'SelfUnion',false);
      % append to running mesh
      SF = [SF;size(SV,1)+SFe];
      J = [J;Je];
      SV = [SV;SVe];
    end
  end
  % Brute force
  [SVAB,SFAB,JAB] = linear_sweep_edges(VA,FA,VB,edges(FB));
  [SVBA,SFBA,JBA] = linear_sweep_edges(VB,FB,VA,edges(FA));
  SV = [SVAB;SVBA];
  SF = [SFAB;size(SVAB,1)+SFBA];
  J =[JAB;size(FA,1)+1+JBA];
  [SV,I] = remove_unreferenced(SV,SF);
  SF = I(SF);
  [SV,~,I] = remove_duplicate_vertices(SV,0);
  SF = I(SF);

  [SV,SF,SJ] = mesh_boolean(SV,SF,[],[],'union');
  J = J(SJ);
end
