%[V,F] = load_mesh('~/Dropbox/models/fandisk.off');
%V = V*axisangle2matrix([1 0 0],pi);
%[V,F] = load_mesh('~/Dropbox/models/bunny.off');
function [tf,te,to,sc,I] = tsurf_cad(F,V,varargin)
  V = V*axisangle2matrix([1 0 0],-pi/2);
  N = normals(V,F);
  BC = barycenter(V,F);

  % sharp edges
  [A,C] = adjacency_dihedral_angle_matrix(V,F);
  %% This is much much slower
  %A(1&A) = abs(A(1&A)-pi)>pi*0.11;
  [AI,AJ,AV] = find(A);
  keep = abs(AV-pi)>(pi*0.11) & ~isnan(AV);
  A = sparse(AI(keep),AJ(keep),1,size(A,1),size(A,2));
  [CI,~,CV] = find(C.*A);
  II = [CI+mod(CV,3)*size(F,1) CI+mod(CV+1,3)*size(F,1)];
  E = F(II);

  % cut mesh at sharp edges to get crisp normals
  [G,I] = cut_edges(F,E);
  W = V(I,:);

  clf;
  hold on;
  blue = [0.2 0.3 0.8];
  tf = tsurf(G,W, ...
    ... % 'FaceVertexCData',repmat(blue,size(W,1),1), ...
    'SpecularStrength',0, ...
    'DiffuseStrength',0.1, ...
    'AmbientStrength',1.0, ...
    'EdgeColor','none','FaceAlpha',0.9,fphong, ...
    varargin{:});
  te = tsurf(E(:,[1 2 2]),V,'EdgeColor',blue*0.75);
  to = tsurf([1 1 1],V,'LineWidth',2,'EdgeColor',blue*0.5);
  view(130,38);
  axis equal;
  l = light('Position',[1 4 5.0],'Style','infinite');
  [h,~,~,g] = add_shadow(tf,l,'Nudge',2e-3,'Fade','local','Color',[0.8 0.8 0.8]);
  % faint amient occlusion
  %AO = ambient_occlusion(W,G,W,per_vertex_normals(W,G),1000);
  %AO = AO*0.27;
  %tf.FaceVertexCData = bsxfun(@times,tf.FaceVertexCData,1-AO);
  %hold off;


  % floor board
  SV = h.Vertices;
  BB = bounding_box(SV(:,1:2));
  M = mean(bounding_box(V(:,1:2)));
  BBmM = bsxfun(@minus,BB,M);
  BB = bsxfun(@plus,sign(BBmM)*max(abs(BBmM(:))),M);
  BB = bsxfun(@plus,bsxfun(@minus,BB,mean(BB))*1.1,mean(BB));
  BB(:,3) = min(V(:,3))-4e-3;
  BB = reshape(permute(BB,[1 3 2]),[2 2 3]);
  % checkboard texture
  nx = 16;
  extent = @(BB) max(max(BB)) - min(min(BB));
  ny = ceil(extent(BB(:,:,2))/ extent(BB(:,:,1))*nx/2)*2;
  ch = repmat(0.99-0.15*xor((mod(repmat(0:8*nx-1,8*ny,1),8*2)>7), ...
    (mod(repmat((0:8*ny-1)',1,8*nx),8*2)>7)),[1 1 3])*0.5 + 0.5;
  hold on;
  sc = surf(BB(:,:,1),BB(:,:,2),BB(:,:,3), ...
    'CData',ch,'FaceColor','texturemap', ...
    'SpecularStrength',0, 'DiffuseStrength',0, 'AmbientStrength',1);
  hold off;

  axis vis3d;
  camproj('persp');
  set(gca,'Visible','off');
  T = get(gca,'tightinset');
  set(gca,'position',[T(1) T(2) 1-T(1)-T(3) 1-T(2)-T(4)]); 

  % Set up rotatation callbacks to hide view-dependent effects during drag
  up = @() ...
    set(to,'Faces', ...
      outline(F((sum(N.*bsxfun(@minus,BC,campos),2)<=0),:))*[1 0 0;0 1 1]) | ...
    set(h,'FaceAlpha',0.5*(g*[campos 1]'<0)) | ...
    set(sc,'FaceAlpha',1.0*(g*[campos 1]'<0));
  up();
  down = @() ...
    (exist('to','var') && ishandle(to) && isempty(set(to,'Faces',[]))) || ...
    isempty(set(rotate3d,'ActionPostCallback',[],'ActionPreCallback',[]));
  set(rotate3d,'ActionPostCallback',@(src,obj) up());
  set(rotate3d,'ActionPreCallback',@(src,obj) down());

  %for t = linspace(0,-360,60)
  %  view(64+t,20);
  %  up();
  %  drawnow;
  %  filename = 'tsurf-cad.gif';
  %  figgif(filename,'nodither');
  %end
end
