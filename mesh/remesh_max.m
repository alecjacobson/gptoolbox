function [VV,FF,J,FJ] = remesh_max(V,F,S)
  % Remesh a triangle mesh so that the index set of the maximum elements of S
  % are piecewise constant over the new triangles.
  %
  % Inputs:
  %   V  #V by dim list of vertex positions
  %   F  #F by 3 list of triangle indices into V
  %   S  #V by #S matrix of scalar values
  % Outputs:
  %   VV  #VV by dim list output vertex positions
  %   FF  #FF by 3 list of output triangle indices
  %   J  #FF list of indices into #S
  %   FJ  #FF list of indices into F
  %

  % V,F as fake input so they are local
  function [VV,FF,J,FJ] = remesh_max_single(S,V,F)
    V = [0 0;1 0;0 1];
    F = [1 2 3];
    FF = F + 3*((1:size(S,2))-1)';
    VV = reshape(permute(repmat(V,[1 1 size(S,2)]),[1 3 2]),[],2);
    FJ = (1:size(S,2))';
    % This is important to ensure/help make single topmost component
    keep = ...
      max([S(FF(:,1)) S(FF(:,2)) S(FF(:,3))],[],2) >= ...
      max(min([S(FF(:,1)) S(FF(:,2)) S(FF(:,3))],[],2)');
    FF = FF(keep,:);
    FJ = FJ(keep);
    [SV,SF,~,SJ] = selfintersect([VV S(:)],FF,'StitchAll',true);
    % This is basically just getting lucky...
    [VV,FF,J] = outer_hull(SV,SF);
    J = FJ(SJ(J));
    %% Elaborate method using extrude
    %[VV,FF] = extrude(VV(:,1:2),FF);
    %% top
    %I =((1:size(FF,1))');
    %I = I<=size(S,2);
    %VV = [VV(:,1:2) [S(:);0*S(:)]];
    %% Could use a Volino type thing (normals+boundary) to cull
    %tic;
    %[VV,FF,J] = mesh_boolean(VV,FF,[],[],'union');
    %toc
    %FF = FF(I(BJ),:);
    %J = J(I(J));
    %[VV,I] = remove_unreferenced(VV,FF);
    %FF=I(FF);
    %% Do I really need to do this just because I don't know SSI? VV(:,3) is already
    %% == max(SS,[],2)
    %B = barycentric_coordinates( ...
    %  VV(:,1:2), ...
    %  V(repmat(1,size(VV,1),1),:), ...
    %  V(repmat(2,size(VV,1),1),:), ...
    %  V(repmat(3,size(VV,1),1),:));
    %SS = B*S;
    %[~,J] = max(SS,[],2);
    %J = FJ(J);
  end

  [~,MI] = max(cat(3,S(F(:,1),:),S(F(:,2),:),S(F(:,3),:)),[],2);
  easy = MI(:,:,1) == MI(:,:,2) & MI(:,:,2) == MI(:,:,3);

  FJ = find(easy);
  J = MI(easy,1,1);
  FF = F(easy,:);
  VV = V;

  not_easy = find(~easy);
  for ne = 1:numel(not_easy)
    progressbar(ne,numel(not_easy));
    f = not_easy(ne);
    Sf = S(F(f,:),:);
    [VVf,FFf,Jf] = remesh_max_single(Sf);
    BC = [1-sum(VVf(:,1:2),2) VVf(:,1) VVf(:,2)];
    VVf = BC * V(F(f,:),:);
    n = size(VV,1);
    VV = [VV;VVf];
    FF = [FF;n+FFf];
    FJ = [FJ;repmat(f,size(FFf,1),1)];
    J = [J;Jf];
    %tsurf(FF,VV,'CData',J);
    %colormap(cbrewer('Set1',size(S,2)));
    %caxis([1 size(S,2)]);
    %view(2);
    %axis equal;
    %drawnow;
  end

end
