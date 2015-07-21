function [G,I] = cut_edges(F,E)
  % CUT_EDGES  Given a list of triangles (F) and a list of edges E to "cut",
  % replicate vertices along cut edges so that they are no longer connected
  % topologically. Vertices are replicated if they no longer share an edge in
  % common. This means that you can cut a single edge inward from the boundary,
  % but you cannot cut just one edge in the interior (you must cut least two)
  %
  % [G,I] = cut_edges(F,E)
  %
  % Inputs:
  %   F  #F by 3 list of triangles indices into V
  %   E  #E by 2 list of edges
  % Outputs:
  %   G  #G by 3 list of triangle indices into W
  %   I  #W by 1 list of indices into V of "birth" vertices
  %   
  % Example:
  %   A = adjacency_matrix(F);
  %   [~,P] = graphshortestpath(A,1);
  %   E = [P{end}(1:end-1);P{end}(2:end)]';
  %   [G,I] = cut_edges(F,E)
  %   W = V(I,:);
  %   tsurf(G,W+10*(massmatrix(W,G)\cotmatrix(W,G)*W))
  %   hold on;
  %   plot_edges(V,E,'r','LineWidth',2);
  %   hold off;
  %

  % List of all "half"-edges: 3*#F by 2
  allE = [F(:,[2 3]); F(:,[3 1]); F(:,[1 2])];
  % Sort each row
  sortallE = sort(allE,2);
  % IC(i) tells us where to find sortallE(i,:) in uE: 
  % so that sortallE(i,:) = uE(IC(i),:)
  [uE,~,IC] = unique(sortallE,'rows');
  % uE2F(e,f) = 1 means face f is adjacent to unique edge e
  uE2F = sparse(IC(:),repmat(1:size(F,1),1,3)',1,size(uE,1),size(F,1));

  m = size(F,1);
  FF = reshape(1:3*m,[],3);
  %VV = V(F,:);

  [~,I] = setdiff(sort(uE,2),sort(E,2),'rows');
  noncut = sparse(I,I,1,size(uE,1),size(uE,1));
  %noncut = speye(size(uE,1));
  uE2EE = sparse(IC(:),(1:3*size(F,1))',1,size(uE,1),3*size(F,1));
  I = reshape(1:3*size(F,1),size(F,1),3);
  VV2EE = sparse(FF(:,[1 2 3 1 2 3]),I(:,[2 3 1 3 1 2]),1,3*size(F,1),3*size(F,1));

  n = max(F(:));
  VV2V = sparse(I,F,1,size(F,1)*3,n);

  A = ...
    (VV2EE * (uE2EE' * noncut * uE2EE) * VV2EE') & ...
    (VV2V * VV2V');

  %[I,J] = find(VV2EE);
  %VV = VV+0.1*(cotmatrix(VV,FF)*VV);
  %allE = [FF(:,[2 3]); FF(:,[3 1]); FF(:,[1 2])];
  %BC = barycenter(VV,allE);
  %uBC = barycenter(V,uE);
  %clf;
  %%tsurf(FF,VV);
  %hold on;
  %[I,J] = find(VV2EE);
  %quiver(VV(I,1),VV(I,2),BC(J,1)-VV(I,1),BC(J,2)-VV(I,2),0)
  %[I,J] = find(uE2EE');
  %quiver(BC(I,1),BC(I,2),uBC(J,1)-BC(I,1),uBC(J,2)-BC(I,2),0,'Color','r')
  %[I,J] = find(VV2EE*(uE2EE'));
  %quiver(VV(I,1),VV(I,2),uBC(J,1)-VV(I,1),uBC(J,2)-VV(I,2),0,'Color','g')
  %%[I,J] = find(VV2V);
  %%quiver(VV(I,1),VV(I,2),V(J,1)-VV(I,1),V(J,2)-VV(I,2),0,'Color','k')
  %[I,J] = find(A);
  %quiver(VV(I,1),VV(I,2),VV(J,1)-VV(I,1),VV(J,2)-VV(I,2),0,'Color',[0.5 0.5 0.5])
  %hold off;
  %full(A)

  %% 1096?
  %A(1094,:)

  I = F(:);
  % dummy vertex data
  VV = zeros(max(FF(:)),1);
  [~,J] = conncomp(A);
  VV(J,:) = VV;
  I(J) = I;
  FF = J(FF);
  [W,IM] = remove_unreferenced(VV,FF);
  I = I(find(IM<=size(W,1)));
  G = IM(FF);

  %tsurf(G,W+4*(massmatrix(W,G)\(cotmatrix(W,G)*W)));
  %tsurf(G,W);
  %spy(uE2EE)
  %tsurf(F,V,'VertexIndices',true);
end
