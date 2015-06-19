function [W,E,U] = skeleton_extraction(V,F)
  % SKELETON_EXTRACTION Compute the "skeleton" of a surface mesh following "Skeleton
  % Extraction by Mesh Contraction" [Au et al. 2008]. The final combinatorially
  % simplification is a bit different than the QSlim-based approach in [Au et
  % al. 2008].
  % 
  % [W,E,U] = skeletonization(V,F)
  % 
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh indices
  % Outputs:
  %   W  #W by 3 list of skeleton vertex positions
  %   E  #E by 2 list of skeleton edges into W
  %   U  #V by 3 list of skeleton vertex positions, before thinning
  %

  delta = 1e-3;
  U = V;
  L = cotmatrix(U,F);
  %tsurf(F,V,'FaceAlpha',0.1,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0.5);
  %hold on;
  %t = tsurf( ...
  %   bsxfun(@plus,size(F,1)*(0:2),(1:size(F,1))'), ...
  %   U(F,:),'FaceColor',[0.7 0.8 1.0],'EdgeColor','k','EdgeAlpha',0.5, ...
  %  'SpecularStrength',0.1,'FaceLighting','phong');
  %tr = tsurf([1 1 1],U,'EdgeColor','r','LineWidth',2);
  %hold off;
  %axis equal;
  %view(-62,10);
  %camproj('persp');
  %light('Position',[-100.0,1.0,-1.0],'Style','infinite');
  %set(gca,'Visible','off');
  %set(gcf,'Color','w');
  %drawnow;

  b = [];

  iter = 0;
  %imwrite(myaa('raw'),sprintf('parasaur-pass-%05d.png',iter));
  %set(t,'EdgeColor','none');
  iter = 1;
  h = avgedge(V,F);

  while true
    L = cotmatrix(U,F);
    M = massmatrix(U,F,'barycentric');
    b = union(find(any(abs(L)>1e+5,2)),b);
    U_prev = U;
    U = min_quad_with_fixed(-(M-delta*L),2*M*U,b,U(b,:));
    %U = U/sqrt(sum(doublearea(U,F))*0.5);
    %U = bsxfun(@minus,U,centroid(U,F));
    %D = max(internalangles(U,F),[],2)>0.9*pi & ...
    %  doublearea(U,F)<1e-5;
    D = all(ismember(F,b),2);
    %set(t,'Vertices',U(F,:));
    %set(tr,'Faces',F(D,:),'Vertices',U);
    %xlabel(sprintf('%d',numel(b)));
    %drawnow;
    if max(max(abs(U-U_prev))) < 1e-2*h
      break;
    end
    %imwrite(myaa('raw'),sprintf('parasaur-pass-%05d.png',iter));
    iter =iter + 1;
  end


  W = U;
  A = adjacency_matrix(F) + speye(size(V,1));
  EC = adjacency_edge_cost_matrix(U,F);
  phi = zeros(size(W,1),1);
  marked = false(size(W,1),1);
  E = [];
  for wi = 1:size(W,1)
    % this is a lot faster than find(A(ui,:))
    if marked(wi)
      continue;
    end
    pI = A(:,wi)>0;
    I = false(size(V,1),1);
    N = [];
    while true
      pN = find(pI);
      if numel(pN) <= 1
        break;
      end
      Wi = W(pN,:);
      [pC,pS,L] = pca(Wi);
      phi(wi) = L(1)/sum(L);
      if phi(wi) < 0.99
        break;
      end
      I = pI;
      N = pN;
      C = pC;
      S = pS;
      %tsurf(F,W);
      %hold on;
      %scatter3(U(I,1),U(I,2),U(I,3),'SizeData',50);
      %hold off;
      %title(sprintf('%g',phi(wi)));
      %axis equal;
      %drawnow;

      %pI = (A*I)>0;
      CCI = sum(A(I,:))';
      ECI = minnz(EC(I,:))';
      ECI = ECI./CCI;
      ECI(I) = 0;
      [~,j] = minnz(ECI);
      pI = I;
      pI(j) = 1;

      if isequal(pI,I)
        break;
      end

    end

    if ~isempty(N)
      [~,O] = sort(S(:,1));
      O = N(O);
      A(I,:) = 0;
      A(:,I) = 0;
      %tsurf(F,W);
      %hold on;
      %scatter3(U(I,1),U(I,2),U(I,3),'SizeData',25);
      %scatter3(U(marked,1),U(marked,2),U(marked,3),'SizeData',50);
      %plot_edges(U,E,'b','LineWidth',3);
      %plot_edges(U,[O(1) O(end)],'r','LineWidth',3);
      %hold off;
      %axis equal;
      %drawnow
      marked(N) = true;
      E = [E;O(1) O(end)];
    end

    %if phi<0.4
    %  continue;
    %end
    %Wi = bsxfun(@plus,bsxfun(@times,S(:,1),C(:,1)'),mean(Wi));
    %d = sqrt(sum(bsxfun(@minus,Wi,Wi(1,:)).^2,2));
    %H = max(d);
    %d = d/H;
    %f = 2.*d.^3 - 3.*d.^2 + 1;
    %f = f/sum(f);
    %W(wi,:) = sum(bsxfun(@times,f,Wi));
  end
    %plot_edges(W,E,'LineWidth',2)

  count = ones(size(W,1),1);
  while true
    [W,IM] = remove_unreferenced(W,E);
    count(IM) = count;
    count = count(1:size(W,1));
    E = IM(E);
    %plot_edges(W,E,'LineWidth',2)
    %axis equal;
    %drawnow;
    A = adjacency_matrix(E)+speye(size(W,1));
    [ncc,C] = conncomp(A);
    if ncc == 1
      break;
    end
    D = pdist2(W,W);
    D(A>0) = inf;
    [D,I] = min(D);
    [~,j] = min(D);
    i = I(j);
    % running average
    W(i,:) = (count(i)*W(i,:)+W(j,:))/(count(i)+1);
    count(i) = count(i) + 1;
    E(E==j) = i;
    E = unique(sort(E,2),'rows');
  end

  %G = F;
  %E = edges(G);
  %
  %A = U(E(:,2),:)-U(E(:,1),:);
  %B = cross(A,U(E(:,1),:),2);
  %UH = [U ones(size(U,1),1)];
  %
  %ne = size(E,1);
  %n = size(U,1);
  %%K = sparse( ...
  %%  repmat(reshape(repmat(1:ne*3,1,3),ne,9),2,1), ...
  %%  bsxfun(@plus,E(:),n*[1 2 0 2 0 1 3 3 3]), ...
  %%  repmat([-A(:,[3 1 2]) A(:,[2 3 1]) -B],2,1), ...
  %%  ne*3,n*4);
  %%Q = K'*K;
  %
  %%K = zeros(3,4,size(E,1));
  %Q = zeros(4,4,size(U,1));
  %%Q = sparse(n*4,n*4);
  %% This can at least be sped up using http://www.alecjacobson.com/weblog/?p=4186
  %for ei = 1:size(E,1)
  %  K = [ ...
  %    0        -A(ei,3)  A(ei,2) -B(ei,1); ...
  %    A(ei,3)  0        -A(ei,1) -B(ei,2); ...
  %    -A(ei,2) A(ei,1)  0        -B(ei,3)];
  %  K2 = K'*K;
  %  Q(:,:,E(ei,1)) = Q(:,:,E(ei,1)) + K2;
  %  Q(:,:,E(ei,2)) = Q(:,:,E(ei,2)) + K2;
  %end
  %% Cost per vertex
  %C = zeros(size(U,1),1);
  %for vi = 1:size(U,1)
  %  C(vi) = UH(vi,:) * (Q(:,:,vi) * (UH(vi,:)'));
  %end
  %
  %
  %A = sparse(repmat((1:ne)',1,2),E,1,ne,n);
  %EL = sqrt(sum(U(E(:,2),:)-U(E(:,1),:),2));
  %% This is a terrible O(n²) implemenation. The energy also doesn't match the [Au
  %% et al. 2008] paper. I don't understand their description or necessarily agree
  %% with the motivation.
  %J = 1:size(U,1);
  %while true
  %  % Edge-vertex incidence matrix
  %  ne = size(E,1);
  %  % Edge lengths summed at vertices
  %  EC = repmat(A*C,2,1) + 0.001*repmat(EL,2,1);
  %  % pick an edge
  %  [~,ei] = min(EC);
  %  e = E(mod(ei-1,ne)+1,:);
  %  if ei>ne
  %    e = fliplr(e);
  %  end
  %  % collapse edge
  %  J(e(2)) = e(1);
  %  J = J(J);
  %  Q(:,:,e(1)) = Q(:,:,e(1)) + Q(:,:,e(2));
  %  C(e(1)) = UH(e(1),:) * (Q(:,:,e(1)) * (UH(e(1),:)'));
  %  G = J(G);
  %  E = J(E);
  %
  %  A(:,e(1)) = A(:,e(1)) + A(:,e(2));
  %  A(:,e(2)) = 0;
  %  ndeg = E(:,1)~=E(:,2);
  %  E = E(ndeg,:);
  %  EL = EL(ndeg,:);
  %  A = A(ndeg,:);
  %  G = G(G(:,1)~=G(:,2)&G(:,2)~=G(:,3)&G(:,3)~=G(:,1),:);
  %  %G = unique(sort(G,2),'rows');
  %  if isempty(G)
  %    break;
  %  end
  %
  %  change = any(E==e(1),2);
  %  Echange = unique(sort(E(change,:),2),'rows');
  %  E = [E(~change,:);Echange];
  %  nc = size(Echange,1);
  %  A = [ ...
  %    A(~change,:); ...
  %    sparse(repmat((1:nc)',1,2),Echange,1,nc,n)];
  %  EL = [EL(~change);sqrt(sum(U(Echange(:,2),:)-U(Echange(:,1),:),2))];
  %
  %  if size(E,1)<50 
  %    tsurf(G,U,'EdgeColor','r');
  %    hold on;
  %    tsurf([E E(:,1)],U,'EdgeColor','k');
  %    hold off;
  %    title(sprintf('%d,%d',size(G,1),size(E,1)));
  %    axis equal;
  %    drawnow
  %  end
  %end
  %M = sparse(J,1:n,diag(massmatrix(V,F)),n,n);
  %M(J,:) = M(J,:)*diag(1./sum(M(J,:),2));
  %U = M*V;
  %[U,IM] = remove_unreferenced(U,E);
  %E = IM(E);
  %plot_edges(U,E,'LineWidth',2);
  %hold on;
  %tsurf(F,V,'CData',J,'FaceAlpha',0.2,'EdgeAlpha',0.2);
  %hold off;
end
