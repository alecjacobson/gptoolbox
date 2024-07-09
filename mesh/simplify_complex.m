function [dE] = simplify_complex(V,E,tol)
  % SIMPLIFY_COMPLEX Simplify a piecewise linear complex of edges by simplifying
  % each connected component between non-manifold edges. Non-manifold vertices
  % are guaranteed to stay in the output.
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   E  #E by 2 list of edge indices into V
  %   tol  dpsimplify tolerance
  % Outputs:
  %   dE  #dE by 2 list of edge indices into V
  % 

  % [dE] = simplify_complex(V,E,tol)
  %
  % for loops, prefer to start at non-manifold vertices
  val = accumarray(E(:),1);
  
  C = manifold_patches(E);
  %tsurf(reshape(1:numel(E),size(E)),V(E,:),'Marker','o','CData',[C;C],'EdgeColor','flat','LineWidth',3);colormap(cbrewer('Set1',max(C)))
  %plot_edges(V,E,'-o','LineWidth',2);
  %txt(barycenter(V,E),num2str(C'),'BackgroundColor',0.9*[1 1 1])
  
  dE = [];
  % implementation is O(n⋅|C|) but "trivial" to factor out if the number of
  % components is ever huge.
  for c = 1:max(C)
    Ec = E(C==c,:);
    [Vc,~,J,Ec] = remove_unreferenced(V,Ec);

    path = edges_to_path(Ec);
    if ~ismember([path(1) path(2)],Ec,'rows')
      path = flip(path);
    end


    is_loop = path(1) == path(end);
    if is_loop
      % there's a loop
      valc = val(J);
      [max_val,i] = max(valc);
      path = path(1:end-1);
      path = path([i:end 1:i-1]);
      assert(path(1) == i);
      path = path([1:end 1]);
    end
    [dVc,dI] = dpsimplify(Vc(path,:),tol);
    dpath = path(dI);
    dEc = [dpath(1:end-1) dpath(2:end)];

    dE = [dE;J(dEc)];
  end
end
